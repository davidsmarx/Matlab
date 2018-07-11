

#define _USE_MATH_DEFINES
#include <math.h>

#include "mex.h"
#include "ProcessRays.h"

#define CM (1.0e-2)
#define MM (1.0e-3)

#define BUFPOS_XYBIN ( 8*sizeof(int) )
#define BUFPOS_X ( 10*sizeof(int) + 2*sizeof(double) )
#define BUFPOS_Y ( 10*sizeof(int) + 3*sizeof(double) )
#define BUFPOS_L ( 10*sizeof(int) + 5*sizeof(double) )
#define BUFPOS_M ( 10*sizeof(int) + 6*sizeof(double) )
#define BUFPOS_INTENSITY ( 10*sizeof(int) + 12*sizeof(double) )
#define BUFPOS_PHASE_AT ( 10*sizeof(int) + 14*sizeof(double) )


void ProcessRays(ProcessRays_Parms *psParms)
{
// for convenience
int Nx = psParms->NPixX;
int Ny = psParms->NPixY;

//
double dx = psParms->dDetectorSizeX/((double)Nx);
double dy = psParms->dDetectorSizeY/((double)Ny);

double *x = new double[Nx]; //dx*(-(Nx/2-1/2):(Nx/2-1/2))';
double *y = new double[Ny]; //dy*(-(Ny/2-1/2):(Ny/2-1/2))';
for ( int ii = 0; ii < Nx; ii++ )
	x[ii] = dx*(-(double)Nx*0.5+0.5 + (double)ii);
for ( int ii = 0; ii < Ny; ii++ )
	y[ii] = dy*(-(double)Ny*0.5+0.5 + (double)ii);

//
//% the zemax rays are x-major
//[Y, X] = meshgrid(y, x);
//% from the zemax manual: For the Detector Rectangle, the pixel numbers start at 1 in the lower left (-x, -y) corner of the rectangle. The
//% pixel numbers increase along the +x direction first, and then move up one row in the +y direction, and start over
//% at the -x edge. For a detector with nx by ny pixels, the pixel numbers are 1 through nx on the bottom row, and
//% then nx+1 through 2nx on the next row above the bottom, until the last pixel (number nx*ny) is reached in the
//% upper right (+x, +y) corner. Another way of stating the pixel numbering is that the x index changes the fastest,
//% and then the y index. check: 
//% % [X(1) Y(1)] is lower left:
//% % [X(N) Y(N)] is lower right:
//% % [X(N+1) Y(N+1)] is left side, up one row (Y)
//% % [X(N*N) Y(N*N)] is upper right:
double *X = new double[Nx*Ny];
double *Y = new double[Nx*Ny];
for ( int ix = 0; ix < Nx; ix++ )
{
	for ( int iy = 0; iy < Ny; iy++ )
	{
		X[ix + iy*Nx] = x[ix];
		Y[ix + iy*Nx] = y[iy];
	}
}


//% calculate phase at each of the four pixels
//% the detector object is rotated 180 about y, but rays.x, rays.y are global
//% coords, but detbin is local. So need to un-rotate local X(bin) back to
//% global coords.
//
//% calculate incoherent irradiance
//numrays = length(rays);
//[IrrMap, AmpMap, FldMap] = deal(zeros(size(X)));
memset(psParms->pdIncIrr, 0, Nx*Ny*sizeof(double));
memset(psParms->pdFieldRe, 0, Nx*Ny*sizeof(double));
memset(psParms->pdFieldIm, 0, Nx*Ny*sizeof(double));
memset(psParms->pdSumAmp, 0, Nx*Ny*sizeof(double));
for ( int iseg = 0; iseg < psParms->Nsegs; iseg++ )
{
	//    thisray = rays(iseg);
	PUINT8 pSeg = psParms->bufferpos_segs[iseg];

	double rayx = (*(double*)(pSeg + BUFPOS_X))*MM;
	double rayy = (*(double*)(pSeg + BUFPOS_Y))*MM;
	double rayl = (*(double*)(pSeg + BUFPOS_L));
	double raym = (*(double*)(pSeg + BUFPOS_M));

	//mexPrintf("rayx = %.3f, rayy = %.3f, rayl = %.3e, raym = %.3e\n", rayx/MM, rayy/MM, rayl, raym);

	//    % looks like xybin is zero-offset, not one-offset
	int detbin = *(int*)(pSeg + BUFPOS_XYBIN); //thisray.xybin+1;

	double xp = X[detbin]; double yp = Y[detbin];
	
	//    % the detector object is rotated 180 about y, rays.x, rays.y are global
	//    % coords, but detbin is local. So need to rotate global rays.x,y about
	//    % y:

	//    xr = -thisray.x;        yr = thisray.y;
	double xr = -rayx; double yr = rayy;

    //% distance from pixel center
    double ddx = xr-xp;
    double ddy = yr-yp;

	//    % double check that ray x,y is closest to bin center given by xybin
	if (fabs(ddx) > 0.5*dx || fabs(ddy) > 0.5*dy)	
		mexErrMsgIdAndTxt("mexReadZRD:ProcessRays", "ray to pixel center check fails");

//    
	int dirx = (int)_copysign(1.0, ddx);
	int diry = (int)_copysign(1.0, ddy);
//    % the other three closest detector pixels
    int detbin_dx = detbin + dirx;
    int detbin_dy = detbin + Nx*diry;
    int detbin_dxdy = detbin + dirx + Nx*diry;

	//    % calculate phase at each of the four pixels
	double k = 2.0*M_PI/psParms->pdWavelength[iseg];
	//        
	//    % ray.phase_at from the zrd file is already calculated at the pixel
	//    % center by zemax. see e-mail exchange with zemax support
	//    Pha = thisray.phase_at + ... 
	//        k.* ( (-X(detbin) - -X(DetBin)).*thisray.l + (Y(detbin) - Y(DetBin)).*thisray.m );
	double pha = *(double*)(pSeg + BUFPOS_PHASE_AT);
	double pha_dx = pha + k*( (-X[detbin] - -X[detbin_dx])*rayl + (Y[detbin] - Y[detbin_dx])*raym );
	double pha_dy = pha + k*( (-X[detbin] - -X[detbin_dy])*rayl + (Y[detbin] - Y[detbin_dy])*raym );
	double pha_dxdy = pha + k*( (-X[detbin] - -X[detbin_dxdy])*rayl + (Y[detbin] - Y[detbin_dxdy])*raym );
	

	//    % calculate field amplitude (=sqrt(intensity)) at each of the four
	//    % pixels	
	//    % bi-linear interpolation to apportion ray value to four pixels
	double ax = fabs(ddx)/dx;
	double ay = fabs(ddy)/dy;

	double rayintensity = *(double*)(pSeg+BUFPOS_INTENSITY);
	double intensity =  (1.0-ax)*(1.0-ay)*rayintensity;
	double intensity_dx = ax*(1.0-ay)*rayintensity;
	double intensity_dy = (1.0-ax)*ay*rayintensity;
	double intensity_dxdy = ax*ay*rayintensity;

	psParms->pdIncIrr[detbin] += intensity;
	psParms->pdIncIrr[detbin_dx] += intensity_dx;
	psParms->pdIncIrr[detbin_dy] += intensity_dy;
	psParms->pdIncIrr[detbin_dxdy] += intensity_dxdy;

	double ampl = sqrt(intensity);
	double ampl_dx = sqrt(intensity_dx);
	double ampl_dy = sqrt(intensity_dy);
	double ampl_dxdy = sqrt(intensity_dxdy);

	// Amp = sqrt(intensity)
	//    Field = Amp .* exp(1i*Pha);
	psParms->pdSumAmp[detbin] += ampl;
	psParms->pdSumAmp[detbin_dx] += ampl_dx;
	psParms->pdSumAmp[detbin_dy] += ampl_dy;
	psParms->pdSumAmp[detbin_dxdy] += ampl_dxdy;

	psParms->pdFieldRe[detbin] += ampl*cos(pha);
	psParms->pdFieldRe[detbin_dx] += ampl_dx*cos(pha_dx);
	psParms->pdFieldRe[detbin_dy] += ampl_dy*cos(pha_dy);
	psParms->pdFieldRe[detbin_dxdy] += ampl_dxdy*cos(pha_dxdy);

	psParms->pdFieldIm[detbin] += ampl*sin(pha);
	psParms->pdFieldIm[detbin_dx] += ampl_dx*sin(pha_dx);
	psParms->pdFieldIm[detbin_dy] += ampl_dy*sin(pha_dy);
	psParms->pdFieldIm[detbin_dxdy] += ampl_dxdy*sin(pha_dxdy);

} // for each ray segment

memcpy(psParms->px, x, Nx*sizeof(double));
memcpy(psParms->py, y, Ny*sizeof(double));

double invdxdy = 1/(dx*dy);
double invsqrtdxdy = sqrt(invdxdy);

// Watts to Watt/m^2 to match units from Zemax
for ( int ii = 0; ii < Nx*Ny; ii++ )
{
	psParms->pdIncIrr[ii] *= invdxdy;
	psParms->pdFieldRe[ii] *= invsqrtdxdy;
	psParms->pdFieldIm[ii] *= invsqrtdxdy;
	psParms->pdSumAmp[ii] *= invsqrtdxdy;
}
// normalize Coherent Field, p. 20 of notebook
for ( int ii = 0; ii < Nx*Ny; ii++ )
{
	if ( psParms->pdSumAmp[ii] > 0 )
	{
		double factor = sqrt(psParms->pdIncIrr[ii]) / psParms->pdSumAmp[ii];
		psParms->pdFieldRe[ii] *= factor;
		psParms->pdFieldIm[ii] *= factor;
	}
}

delete [] x;
delete [] y;
delete [] X;
delete [] Y;

} // ProcessRays