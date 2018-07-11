#ifndef PROCESSRAYS_H
#define PROCESSRAYS_H

#include <Windows.h>

typedef struct
{
	// inputs
	double dDetectorSizeX, dDetectorSizeY;
	int NPixX, NPixY;
	int Nsegs;
	PUINT8 *bufferpos_segs;
	double *pdWavelength; // [segsfoundcnt]

	// outputs:
	// allocated to NPixX x NPixY arrays in the calling routine
	double *pdIncIrr, *pdFieldRe, *pdFieldIm, *pdSumAmp;
	double *px, *py;

} ProcessRays_Parms;

void ProcessRays(ProcessRays_Parms *psParms);


#endif