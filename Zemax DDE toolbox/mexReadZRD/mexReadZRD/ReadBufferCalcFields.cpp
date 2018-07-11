
#include <windows.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "mex.h"

#include "ReadBufferCalcFields.h"

using namespace std;

#define CM (1.0e-2)
#define MM (1.0e-3)
#define UM (1.0e-6)

#define BUFPOS_FIRSTSEGMENT (sizeof(int))
#define SEGPOS_HIT_OBJECT (2*sizeof(int))
#define SEGPOS_XYBIN ( 8*sizeof(int) )
#define SEGPOS_WAVELENGTH (10*sizeof(long)+8*sizeof(double))
#define SEGPOS_X ( 10*sizeof(int) + 2*sizeof(double) )
#define SEGPOS_Y ( 10*sizeof(int) + 3*sizeof(double) )
#define SEGPOS_L ( 10*sizeof(int) + 5*sizeof(double) )
#define SEGPOS_M ( 10*sizeof(int) + 6*sizeof(double) )
#define SEGPOS_INTENSITY ( 10*sizeof(int) + 12*sizeof(double) )
#define SEGPOS_PHASE_AT ( 10*sizeof(int) + 14*sizeof(double) )

class CDetectorFields
{
	int m_Nx, m_Ny;
	double m_dx, m_dy;
	double *m_x, *m_y, *m_X, *m_Y;

	double *m_pdIncIrr, *m_pdFieldRe, *m_pdFieldIm, *m_pdSumAmp;

public:

	CDetectorFields(ReadBufferCalcFields_Parms *psParms)
	{
		m_pdIncIrr = psParms->pdIncIrr;
		m_pdFieldRe = psParms->pdFieldRe;
		m_pdFieldIm = psParms->pdFieldIm;
		m_pdSumAmp = psParms->pdSumAmp;

		// initialize detector array pixels and fields
		// for convenience
		m_Nx = psParms->NPixX;
		m_Ny = psParms->NPixY;

		//
		m_dx = psParms->dDetectorSizeX/((double)m_Nx);
		m_dy = psParms->dDetectorSizeY/((double)m_Ny);

		m_x = new double[m_Nx]; //dx*(-(m_Nx/2-1/2):(m_Nx/2-1/2))';
		m_y = new double[m_Ny]; //dy*(-(m_Ny/2-1/2):(m_Ny/2-1/2))';
		m_X = new double[m_Nx*m_Ny];
		m_Y = new double[m_Nx*m_Ny];


		for ( int ii = 0; ii < m_Nx; ii++ )
			m_x[ii] = m_dx*(-(double)m_Nx*0.5+0.5 + (double)ii);
		for ( int ii = 0; ii < m_Ny; ii++ )
			m_y[ii] = m_dy*(-(double)m_Ny*0.5+0.5 + (double)ii);

		//
		//% the zemax rays are x-major
		//[Y, X] = meshgrid(y, x);
		//% from the zemax manual: For the Detector Rectangle, the pixel numbers start at 1 in the lower left (-x, -y) corner of the rectangle. The
		//% pixel numbers increase along the +x direction first, and then move up one row in the +y direction, and start over
		//% at the -x edge. For a detector with m_Nx by m_Ny pixels, the pixel numbers are 1 through m_Nx on the bottom row, and
		//% then m_Nx+1 through 2m_Nx on the next row above the bottom, until the last pixel (number m_Nx*m_Ny) is reached in the
		//% upper right (+x, +y) corner. Another way of stating the pixel numbering is that the x index changes the fastest,
		//% and then the y index. check: 
		//% % [X(1) Y(1)] is lower left:
		//% % [X(N) Y(N)] is lower right:
		//% % [X(N+1) Y(N+1)] is left side, up one row (Y)
		//% % [X(N*N) Y(N*N)] is upper right:
		for ( int ix = 0; ix < m_Nx; ix++ )
		{
			for ( int iy = 0; iy < m_Ny; iy++ )
			{
				m_X[ix + iy*m_Nx] = m_x[ix];
				m_Y[ix + iy*m_Nx] = m_y[iy];
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
		// these arrays are allocated in mexReadZRD()
		memset(m_pdIncIrr, 0, m_Nx*m_Ny*sizeof(double));
		memset(m_pdFieldRe, 0, m_Nx*m_Ny*sizeof(double));
		memset(m_pdFieldIm, 0, m_Nx*m_Ny*sizeof(double));
		memset(m_pdSumAmp, 0, m_Nx*m_Ny*sizeof(double));
	}  // INitializeDetector

	~CDetectorFields()
	{
		delete [] m_x;
		delete [] m_y;
		delete [] m_X;
		delete [] m_Y;
	}

	void UpdateDetectorFields(double dWavelength, PUINT8 pSeg)
	{
		////    thisray = rays(iseg);
		//PUINT8 pSeg = psParms->bufferpos_segs[iseg];

		double rayx = (*(double*)(pSeg + SEGPOS_X))*MM;
		double rayy = (*(double*)(pSeg + SEGPOS_Y))*MM;
		double rayl = (*(double*)(pSeg + SEGPOS_L));
		double raym = (*(double*)(pSeg + SEGPOS_M));

		//mexPrintf("rayx = %.3f, rayy = %.3f, rayl = %.3e, raym = %.3e\n", rayx/MM, rayy/MM, rayl, raym);

		//    % looks like xybin is zero-offset, not one-offset
		int detbin = *(int*)(pSeg + SEGPOS_XYBIN); //thisray.xybin+1;

		double xp = m_X[detbin]; double yp = m_Y[detbin];
	
		//    % the detector object is rotated 180 about y, rays.x, rays.y are global
		//    % coords, but detbin is local. So need to rotate global rays.x,y about
		//    % y:

		//    xr = -thisray.x;        yr = thisray.y;
		double xr = -rayx; double yr = rayy;

		//% distance from pixel center
		double ddx = xr-xp;
		double ddy = yr-yp;

		////    % double check that ray x,y is closest to bin center given by xybin
		//if (fabs(ddx) > 0.5*m_dx || fabs(ddy) > 0.5*m_dy)	
		//	mexErrMsgIdAndTxt("mexReadZRD:ProcessRays", "ray to pixel center check fails");

	//    
		int dirx = (int)_copysign(1.0, ddx);
		int diry = (int)_copysign(1.0, ddy);
	//    % the other three closest detector pixels
		int detbin_dx = detbin + dirx;
		int detbin_dy = detbin + m_Nx*diry;
		int detbin_dxdy = detbin + dirx + m_Nx*diry;

		//    % calculate phase at each of the four pixels
		double k = 2.0*M_PI/dWavelength;
		//        
		//    % ray.phase_at from the zrd file is already calculated at the pixel
		//    % center by zemax. see e-mail exchange with zemax support
		//    Pha = thisray.phase_at + ... 
		//        k.* ( (-X(detbin) - -X(DetBin)).*thisray.l + (Y(detbin) - Y(DetBin)).*thisray.m );
		double pha = *(double*)(pSeg + SEGPOS_PHASE_AT);
		double pha_dx = pha + k*( (-m_X[detbin] - -m_X[detbin_dx])*rayl + (m_Y[detbin] - m_Y[detbin_dx])*raym );
		double pha_dy = pha + k*( (-m_X[detbin] - -m_X[detbin_dy])*rayl + (m_Y[detbin] - m_Y[detbin_dy])*raym );
		double pha_dxdy = pha + k*( (-m_X[detbin] - -m_X[detbin_dxdy])*rayl + (m_Y[detbin] - m_Y[detbin_dxdy])*raym );
	

		//    % calculate field amplitude (=sqrt(intensity)) at each of the four
		//    % pixels	
		//    % bi-linear interpolation to apportion ray value to four pixels
		double ax = fabs(ddx)/m_dx;
		double ay = fabs(ddy)/m_dy;

		double rayintensity = *(double*)(pSeg+SEGPOS_INTENSITY);
		double intensity =  (1.0-ax)*(1.0-ay)*rayintensity;
		double intensity_dx = ax*(1.0-ay)*rayintensity;
		double intensity_dy = (1.0-ax)*ay*rayintensity;
		double intensity_dxdy = ax*ay*rayintensity;

		m_pdIncIrr[detbin] += intensity;
		m_pdIncIrr[detbin_dx] += intensity_dx;
		m_pdIncIrr[detbin_dy] += intensity_dy;
		m_pdIncIrr[detbin_dxdy] += intensity_dxdy;

		double ampl = sqrt(intensity);
		double ampl_dx = sqrt(intensity_dx);
		double ampl_dy = sqrt(intensity_dy);
		double ampl_dxdy = sqrt(intensity_dxdy);

		// Amp = sqrt(intensity)
		//    Field = Amp .* exp(1i*Pha);
		m_pdSumAmp[detbin] += ampl;
		m_pdSumAmp[detbin_dx] += ampl_dx;
		m_pdSumAmp[detbin_dy] += ampl_dy;
		m_pdSumAmp[detbin_dxdy] += ampl_dxdy;

		m_pdFieldRe[detbin] += ampl*cos(pha);
		m_pdFieldRe[detbin_dx] += ampl_dx*cos(pha_dx);
		m_pdFieldRe[detbin_dy] += ampl_dy*cos(pha_dy);
		m_pdFieldRe[detbin_dxdy] += ampl_dxdy*cos(pha_dxdy);

		m_pdFieldIm[detbin] += ampl*sin(pha);
		m_pdFieldIm[detbin_dx] += ampl_dx*sin(pha_dx);
		m_pdFieldIm[detbin_dy] += ampl_dy*sin(pha_dy);
		m_pdFieldIm[detbin_dxdy] += ampl_dxdy*sin(pha_dxdy);

	} // UpdateDetectorFields

	void NormalizeDetectorFields(void)
	{
		double invdxdy = 1/(m_dx*m_dy);
		double invsqrtdxdy = sqrt(invdxdy);

		// Watts to Watt/m^2 to match units from Zemax
		for ( int ii = 0; ii < m_Nx*m_Ny; ii++ )
		{
			m_pdIncIrr[ii] *= invdxdy;
			m_pdFieldRe[ii] *= invsqrtdxdy;
			m_pdFieldIm[ii] *= invsqrtdxdy;
			m_pdSumAmp[ii] *= invsqrtdxdy;
		}
		// normalize Coherent Field, p. 20 of notebook
		for ( int ii = 0; ii < m_Nx*m_Ny; ii++ )
		{
			if ( m_pdSumAmp[ii] > 0 )
			{
				double factor = sqrt(m_pdIncIrr[ii]) / m_pdSumAmp[ii];
				m_pdFieldRe[ii] *= factor;
				m_pdFieldIm[ii] *= factor;
			}
		}

	} // NormalizeDetectorFields

}; // class DetectorFields

//void ReadBufferCalcFields(ReadBufferCalcFields_Parms *psParms)
DWORD WINAPI ReadBufferCalcFields (LPVOID lpParam)
{
	ReadBufferCalcFields_Parms* psParms = (ReadBufferCalcFields_Parms*) lpParam;

	UINT8 *buffer = psParms->buffer; // for convenience

    ULONG32 Raypos = 0; // points to a ray (int32 = # of segments to follow) in the buffer
    ULONG32 Segpos = 0; // points to a ray segment in the buffer
    
	// array of pointers to the buffer to record list of ray segments to be returned
	ULONG32 lSizeofRaypath = sizeof(RAYPATH); // the length of a ray segment record in the binary buffer

	CDetectorFields DetectorFields(psParms);
   
	//long version = *(long*)(buffer);
	//long maxnumsegs = *(long*)(buffer+4);
    Raypos = 8; // first ray
	long sourcecnt = 0;
	long segsfoundcnt = 0;
	while ( (Raypos < psParms->buffersize) ) //&& (sourcecnt < lMaxRays) )
	{
		// # of segments for this ray
        long numsegsthisray = *(long*)(buffer+Raypos);

		Segpos = Raypos + BUFPOS_FIRSTSEGMENT; // move to the ray segment
        
		//Segpos += 8; // move to the "hit_object" member of the struct, which is the source
		long src_object = *(long*)(buffer+Segpos+SEGPOS_HIT_OBJECT);

		// test that this ray has enough segments to get to the first segment that might hit the detector
		// select only rays that start from SOURCE_OBJECT, since this is the first segment of the ray, its hit_object is the source
		if ( src_object == psParms->iSourceObj && numsegsthisray >= psParms->lFirstSegment )
		{
			// wavelength is stored in segment 0
			double dWavetmp = (*(double*)(buffer+Segpos+SEGPOS_WAVELENGTH))*UM;

			// check if this ray hits the detector object at a later segment
			Segpos += psParms->lFirstSegment*lSizeofRaypath;  // advance to the first segment after the source
			long segcnt = 1; // count the number of ray segments until it hits the detecotr
			for ( ; Segpos < Raypos + numsegsthisray*lSizeofRaypath; Segpos += lSizeofRaypath )
			{

				// if this ray segment hits the detector
				if ( (*(long*)(buffer+Segpos+SEGPOS_HIT_OBJECT)) == psParms->lDetectorObject )
				{

					DetectorFields.UpdateDetectorFields(dWavetmp, buffer+Segpos);
	
					segsfoundcnt++;

					//break; // done, go to the next ray
				}

				segcnt++;
			} // for each segment in this ray
		} // if ray is from requested source

        // advance the Raypos to the next launched ray
		Raypos += BUFPOS_FIRSTSEGMENT + numsegsthisray*lSizeofRaypath;
		
		// increment the ray counter
		sourcecnt ++;
	} // for each launched ray

	//mexPrintf("number of source rays = %d, number of segments found = %d, last source object = %d\n", sourcecnt, segsfoundcnt, src_object[sourcecnt-1]);

	// done, psParms->pdIncIrr, etc. contain the calculated fields for this source and detector
	DetectorFields.NormalizeDetectorFields();

	// other outputs:
	*(psParms->pNumSegsFound) = segsfoundcnt;

	return 1;

} // ReadBufferCalcFields
