#ifndef READBUFFERCALCFIELDS_H
#define READBUFFERCALCFIELDS_H

#include <Windows.h>

typedef struct
{
	// inputs
	INT32 iSourceObj;
	double dDetectorSizeX, dDetectorSizeY;
	int NPixX, NPixY;
	size_t buffersize;
	UINT8 *buffer;
	long lFirstSegment; // skip ray segments to lFirstSegment when searching for the segment that hits the detector
	long lDetectorObject; // 

	// outputs:
	int *pNumSegsFound;
	// allocated to NPixX x NPixY arrays in the calling routine
	double *pdIncIrr, *pdFieldRe, *pdFieldIm, *pdSumAmp;
	double *px, *py;

} ReadBufferCalcFields_Parms;

typedef struct
{
	unsigned int status;
	int level;
	int hit_object;
	int hit_face;
	int unused;
	int in_object;
	int parent;
	int storage;
	int xybin, lmbin; // 10 int32
	double index, starting_phase;
	double x, y, z;
	double l, m, n;
	double nx, ny, nz;
	double path_to, intensity;
	double phase_of, phase_at;
	double exr, exi, eyr, eyi, ezr, ezi; // 21 doubles
} RAYPATH;

//void ReadBufferCalcFields(ReadBufferCalcFields_Parms *psParms);
DWORD WINAPI ReadBufferCalcFields (LPVOID lpParam);

#endif