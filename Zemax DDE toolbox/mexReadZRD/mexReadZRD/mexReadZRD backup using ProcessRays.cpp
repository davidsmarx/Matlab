/*=================================================================
 * mxcreatestructarray.c
 *
 * mxcreatestructarray illustrates how to create a MATLAB structure
 * from a corresponding C structure.  It creates a 1-by-4 structure mxArray,
 * which contains two fields, name and phone number where name is store as a
 * string and phone number is stored as a double.  The structure that is
 * passed back to MATLAB could be used as input to the phonebook.c example
 * in $MATLAB/extern/examples/refbook.
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2011 The MathWorks, Inc.
 * All rights reserved.
 *=================================================================*/

#include <windows.h>
#include <math.h>
#include "mex.h"

#include "ProcessRays.h"

using namespace std;

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

void ReadSegment(UINT8 *bufferSeg, mxArray **pmxIntData, mxArray **pmxDoubleData)
{
    // bufferSeg points to the start of a segment (status field)
    // allcoate mxArray to hold all the int's
    // allocate mxArray to hold all the doubles  in this segment
    
    // for debug:
    mexPrintf("Segment:\nStatus = %u\nhit_object = %d\n",
            *(unsigned int*)(bufferSeg)
            ,*(int*)(bufferSeg+8));
    
    long *pintData; pintData = (long*)mxCalloc(10, sizeof(long));
    double *pdoubleData = (double*)mxCalloc(21, sizeof(double));
    memset(pintData, 0, 10*sizeof(long));
    memset(pdoubleData, 0, 21*sizeof(double));

    memcpy(pintData, bufferSeg, 10*sizeof(long));
    memcpy(pdoubleData, bufferSeg+10*sizeof(long), 21*sizeof(double));
    
    /* Create a 0-by-0 mxArray; you will allocate the memory dynamically */
    *pmxIntData = mxCreateNumericMatrix(0, 0, mxUINT32_CLASS, mxREAL);
    *pmxDoubleData = mxCreateDoubleMatrix(0, 0, mxREAL);
    
    /* Point mxArray to dynamicData */
    mxSetData(*pmxIntData, pintData);
    mxSetM(*pmxIntData, 10);
    mxSetN(*pmxIntData, 1);
    mxSetData(*pmxDoubleData, pdoubleData);
    mxSetM(*pmxDoubleData, 21); 
    mxSetN(*pmxDoubleData, 1);

} // ReadSegment

void CopyArrayToOutput(long lN, double *pdArray, mxArray **pmxDoubleArray)
{
	/* Create an m-by-n mxArray; you will copy existing data into it */
    *pmxDoubleArray = mxCreateDoubleMatrix(lN, 1, mxREAL);
    double *pointer = mxGetPr(*pmxDoubleArray);

	memcpy(pointer, pdArray, lN*sizeof(double));
}

void CopyArrayToOutput(long lN, long *plArray, mxArray **pmxDataArray)
{
    *pmxDataArray = mxCreateNumericMatrix(lN, 1, mxINT32_CLASS, mxREAL);
	long *pointer = (long*)mxGetData(*pmxDataArray);

	memcpy(pointer, plArray, lN*sizeof(long));
}

void Segment2mxStruct(long lNumSegs, PUINT8 bufferSeg[], double pdWavelength[], mxArray **pmxRayStruct)
{
	// field names are in exact order of RAYPATH struct as defined in the zemax manual so that we can
	// copy by index 1-to-1 from binary ray segment to the matlab struct. Then "wavelength" is added on the end.
	const char *field_names[] = { "status"
	,"level"
	,"hit_object"
	,"hit_face"
	,"unused"
	,"in_object"
	,"parent"
	,"storage"
	,"xybin", "lmbin"
	,"index", "starting_phase"
	,"x", "y", "z" // units (mm)
	,"l", "m", "n"
	,"nx", "ny", "nz"
	,"path_to", "intensity"
	,"phase_of", "phase_at"
	,"exr", "exi", "eyr", "eyi", "ezr", "ezi"
	,"wavelength"};
	const long lNumFields = (sizeof(field_names)/sizeof(*field_names));
	const long lNumInts = 10;
	const long lNumDbls = 21;
	const long lPosDbls = lNumInts*sizeof(long); // how many bytes into the record for the double data
	const double dMM = 1.0e-3; // MM
	const double dUM = 1.0e-6; // UM

	//mexPrintf("# of field names = %d\n", lNumFields);

	// create an array of struct lNumSegs x 1
	mwSize dims[2] = {lNumSegs, 1};
	*pmxRayStruct = mxCreateStructArray(2, dims, lNumFields, field_names);

	 /* Since we just
       created the structure and the field number indices are zero
       based, status will always be field# 0 and level will always
       be field# 1, etc. */

	for ( int iseg = 0; iseg < lNumSegs; iseg++ )
	{
		// read in the ints and then the doubles for this ray segment
		for ( int ii = 0; ii < lNumInts; ii++ )
		{
			long ltmp = *(long*)(bufferSeg[iseg]+ii*sizeof(long));
			mxArray* mxDbl = mxCreateDoubleScalar((double)ltmp); // allocates memory
			//*mxGetPr(mxDbl) = (double)ltmp;
			mxSetFieldByNumber(*pmxRayStruct, iseg, ii, mxDbl); /// assigns this field to point to the memory
		}

		// read in the double value data
		int ii = 0;
		for ( ii = 0; ii < 2; ii++ )
		{
			mxArray* mxDbl = mxCreateDoubleScalar( *(double*)(bufferSeg[iseg] + lPosDbls + ii*sizeof(double)) ); // allocates memory
			mxSetFieldByNumber(*pmxRayStruct, iseg, lNumInts+ii, mxDbl);
		}
		for ( ; ii < 5; ii++ ) 		// mks units, x, y, z are in mm
		{
			mxArray* mxDbl = mxCreateDoubleScalar(dMM * (*(double*)(bufferSeg[iseg] + lPosDbls + ii*sizeof(double))) ); // allocates memory
			mxSetFieldByNumber(*pmxRayStruct, iseg, lNumInts+ii, mxDbl);
		}
		for ( ; ii < lNumDbls; ii++ ) // remainder of the double value data
		{
			mxArray* mxDbl = mxCreateDoubleScalar( *(double*)(bufferSeg[iseg] + lPosDbls + ii*sizeof(double)) ); // allocates memory
			mxSetFieldByNumber(*pmxRayStruct, iseg, lNumInts+ii, mxDbl);
		}
		// add wavelength to the end, wavelength units = um
		mxSetFieldByNumber(*pmxRayStruct, iseg, lNumInts+ii, mxCreateDoubleScalar(pdWavelength[iseg])); // already mks


	} // for each segment

} // Segment2mxStruct

void CopyTranspose(double *pIn, double *pOut, int Nx, int Ny)
{
	for ( int ix = 0; ix < Nx; ix++ )
		for ( int iy = 0; iy < Ny; iy++ )
			pOut[iy + Ny*ix] = pIn[ix + Nx*iy];
} // CopyTranspose

void OutputMapsToStruct(ProcessRays_Parms *pSMaps,  mxArray **pmxMapStruct)
{	
	const char *field_names[] = { "IncIrr", "Field", "SumAmp", "x", "y" };
	const long lNumFields = (sizeof(field_names)/sizeof(*field_names));

	// create struct 1x1 array
	mwSize dims[2] = {1, 1};
	*pmxMapStruct = mxCreateStructArray(2, dims, lNumFields, field_names);
	/* Since we just
       created the structure and the field number indices are zero
       based, , etc. */
	
	mxArray *IncIrr = mxCreateDoubleMatrix(pSMaps->NPixX, pSMaps->NPixY, mxREAL);
	double *pointer = mxGetPr(IncIrr);
   // memcpy(pointer, pSMaps->pdIncIrr, pSMaps->NPixX * pSMaps->NPixY * sizeof(double));
	CopyTranspose(pSMaps->pdIncIrr, pointer, pSMaps->NPixX, pSMaps->NPixY);
	mxSetFieldByNumber(*pmxMapStruct, 0, 0, IncIrr);

	mxArray *CohField = mxCreateDoubleMatrix(pSMaps->NPixX, pSMaps->NPixY, mxCOMPLEX);
	double *pRe = mxGetPr(CohField);
	double *pIm = mxGetPi(CohField);
	//memcpy(pRe, pSMaps->pdFieldRe,  pSMaps->NPixX * pSMaps->NPixY * sizeof(double));
	CopyTranspose(pSMaps->pdFieldRe, pRe, pSMaps->NPixX, pSMaps->NPixY);
	//memcpy(pIm, pSMaps->pdFieldIm,  pSMaps->NPixX * pSMaps->NPixY * sizeof(double));
	CopyTranspose(pSMaps->pdFieldIm, pIm, pSMaps->NPixX, pSMaps->NPixY);
	mxSetFieldByNumber(*pmxMapStruct, 0, 1, CohField);

	mxArray *SumAmp = mxCreateDoubleMatrix(pSMaps->NPixX, pSMaps->NPixY, mxREAL);
	double *pSA = mxGetPr(SumAmp);
	CopyTranspose(pSMaps->pdSumAmp, pSA, pSMaps->NPixX, pSMaps->NPixY);
	mxSetFieldByNumber(*pmxMapStruct, 0, 2, SumAmp);

	mxArray *x = mxCreateDoubleMatrix(pSMaps->NPixX, 1, mxREAL );
	double *px = mxGetPr(x);
	memcpy(px, pSMaps->px, pSMaps->NPixX * sizeof(double));
	mxSetFieldByNumber(*pmxMapStruct, 0, 3, x);

	mxArray *y = mxCreateDoubleMatrix(pSMaps->NPixY, 1, mxREAL );
	double *py = mxGetPr(y);
	memcpy(py, pSMaps->py, pSMaps->NPixY * sizeof(double));
	mxSetFieldByNumber(*pmxMapStruct, 0, 4, y);



} // OutputMapsToStruct

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
    
    /* Check for proper number of input and  output arguments */    
    if (nrhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:mxcreatestructarray:maxrhs",
                "usage: mexReadZRD(filename, sOptions)");
    } 

	bool bReturnRays = false;
    switch (nlhs )
	{
	case 1:
		bReturnRays = false;
		break;
	case 2:
		bReturnRays = true;
		break;
	default:
        mexErrMsgIdAndTxt( "MATLAB:mxcreatestructarray:maxlhs",
                "Too many output arguments.");
    }

	// input arguments
    size_t buffersize = mxGetNumberOfElements(prhs[0]);
   // mexPrintf("number of bytes = %u\n",buffersize);
    
    UINT8 *buffer = (UINT8*)mxGetData(prhs[0]);
    
	// input parameters sOptions
    if ( !mxIsStruct(prhs[1]) )
	{
		mexErrMsgIdAndTxt( "MATLAB:mxcreatestructarray:maxrhs",
			"sOptions must be a struct");
	}

	//////////////////////// Validate Source Ojbect List and build Index table
	struct
	{
		size_t N;
		INT32* pList;
	} sSourceList;

	int iFieldSourceObject = mxGetFieldNumber(prhs[1], "SourceObject");
	if ( iFieldSourceObject == -1 ) mexErrMsgIdAndTxt( "mexReadZRD:SourceObject", "No Source Object!");

	mxArray *mxtmp = mxGetField(prhs[1], 0, "SourceObject");
	if ( mxGetClassID(mxtmp) != mxINT32_CLASS )
	{
		mexErrMsgIdAndTxt( "mexReadZRD:SourceObject", "SourceObject must be int32" );
	}
	sSourceList.N = mxGetNumberOfElements(mxtmp);
	if ( sSourceList.N <= 0 ) mexErrMsgIdAndTxt( "mexReadZRD:SourceObject", "SourceObject cannot be empty" );
	sSourceList.pList = (INT32*)mxGetData(mxtmp);

	// find maximum source object # so we can build an index table
	INT32 lMaxSObjNum = 0;
	for ( int ii = 0; ii < sSourceList.N; ii++ )
		if ( sSourceList.pList[ii] > lMaxSObjNum ) lMaxSObjNum = sSourceList.pList[ii];

	int *pSIT = (int*)mxCalloc(lMaxSObjNum, sizeof(int)); // Source Index Table
	for ( int ii = 0; ii < sSourceList.N; ii++ )
		pSIT[sSourceList.pList[ii]] = ii; // ii wil be the index into arrays for storing data for source pList[ii]

	//mexPrintf("Number of SourceObjects is %d, first source object = %d, Max Source Obj Num = %d\n", sSourceList.N, sSourceList.pList[0], lMaxSObjNum);


	////// validate Detector Object
	long lDetectorObject = 46; // default value
	int iFieldDetectorObject = mxGetFieldNumber(prhs[1], "DectectorObject");
	if ( iFieldDetectorObject != -1 )
	{
		mxArray *mxtmp = mxGetField(prhs[1], 0, "DectectorObject");
		if ( mxIsDouble(mxtmp) )		
			lDetectorObject = (long)(*mxGetPr(mxtmp));
	}

	/////// Validate First Segment field
	long lFirstSegment = 1; // default value, first segment to begin search for detector
	int iFieldFirstSegment = mxGetFieldNumber(prhs[1], "FirstSegmentToSearch");
	if ( iFieldFirstSegment != -1 )
	{
		mxArray *mxtmp = mxGetField(prhs[1], 0, "FirstSegmentToSearch");
		if ( mxIsDouble(mxtmp) )		
			lFirstSegment = (long)(*mxGetPr(mxtmp));
	}

	//mexPrintf("input source object = %d\n",lSourceObject);

	////////// validate Detector Size x,y
	double *pdDetectorSize;
	int iFieldDetSize = mxGetFieldNumber(prhs[1],"DetectorSize");
	if ( iFieldDetSize == -1 ) mexErrMsgIdAndTxt( "mexReadZRD:DetectorSize", "DetectorSize missing");
	mxtmp = mxGetField(prhs[1], 0, "DetectorSize");
	if ( !mxIsDouble(mxtmp) || (mxGetNumberOfElements(mxtmp) != 2) ) mexErrMsgIdAndTxt( "mexReadZRD:DetectorSize", "DetectorSize is not double of length 2");
	pdDetectorSize = mxGetPr(mxtmp);

	////////// validate number of pixels x,y
	double *pdNPix;
	int iFieldDetPix = mxGetFieldNumber(prhs[1],"DetectorPixels");
	if ( iFieldDetPix == -1 ) mexErrMsgIdAndTxt( "mexReadZRD:DetectorPixels", "DetectorPixels missing");
	mxtmp = mxGetField(prhs[1], 0, "DetectorPixels");
	if ( !mxIsDouble(mxtmp) || (mxGetNumberOfElements(mxtmp) != 2) ) mexErrMsgIdAndTxt( "mexReadZRD:DetectorPixels", "DetectorPixels is not double of length 2");
	pdNPix = mxGetPr(mxtmp);



    ULONG32 Raypos = 0; // points to a ray (int32 = # of segments to follow) in the buffer
    ULONG32 Segpos = 0; // points to a ray segment in the buffer
    
	// parse the binary zrd buffer
	long version = *(long*)(buffer);
	long maxnumsegs = *(long*)(buffer+4);
    
	//mexPrintf("version = %d, max num of segs = %d\n", version, maxnumsegs );

	// allocate an array and count number of starting rays and their "hit_object"
	const long lMaxRays = 1000000;

	// arrays to count statistics
	long *numsegsperray = (long*)mxCalloc(lMaxRays, sizeof(long)); // store the number of segs for each source ray
	long *src_object = (long*)mxCalloc(lMaxRays, sizeof(long)); // store the hit_object (source obj) for each source ray
	double *pdWavelength = (double*)mxCalloc(lMaxRays, sizeof(long)); // store wavelength for each ray launched

	// array of pointers to the buffer to record list of ray segments to be returned
	ULONG32 lSizeofRaypath = sizeof(RAYPATH); // the length of a ray segment record in the binary buffer
	PUINT8 *bufferpos_rayseg = (PUINT8*)mxCalloc(lMaxRays, sizeof(PUINT8));

	// initialize all arrays to zero
	long lhitobj = 0;
	memset(numsegsperray, 0, lMaxRays*sizeof(long));
	memset(src_object, 0, lMaxRays*sizeof(long));
	memset(bufferpos_rayseg, 0, lMaxRays*sizeof(PUINT8));
   

    Raypos = 8; // first ray
	long sourcecnt = 0;
	long segsfoundcnt = 0;
	while ( (Raypos < buffersize) && (sourcecnt < lMaxRays) )
	{
		// # of segments for this ray
        numsegsperray[sourcecnt] = *(long*)(buffer+Raypos);

		Segpos = Raypos + 4; // move to the ray segment
        
		//Segpos += 8; // move to the "hit_object" member of the struct, which is the source
		src_object[sourcecnt] = *(long*)(buffer+Segpos+8);

		// select only rays that start from SOURCE_OBJECT, since this is the first segment of the ray, its hit_object is the source
		if ( src_object[sourcecnt] == sSourceList.pList[0] )
		{
			// wavelength is stored in segment 0
			double dWavetmp = *(double*)(buffer+Segpos+10*sizeof(long)+8*sizeof(double));

			// check if this ray hits the detector object at a later segment
			Segpos += lFirstSegment*lSizeofRaypath;  // advance to the first segment after the source
			long segcnt = 1; // count the number of ray segments until it hits the detecotr
			for ( ; Segpos < Raypos + numsegsperray[sourcecnt]*lSizeofRaypath; Segpos += lSizeofRaypath )
			{

				lhitobj = *(long*)(buffer+Segpos+8);
				if ( lhitobj == lDetectorObject )
				{

					// want to use this segment
					bufferpos_rayseg[segsfoundcnt] = buffer+Segpos;
					pdWavelength[segsfoundcnt] = dWavetmp*1.0e-6; // um
						
	
					segsfoundcnt++;

					//break; // done, go to the next ray
				}

				segcnt++;
			} // for each segment in this ray
		} // if ray is from requested source

        // move the remainder of this segment + the rest of the segments for this ray
		Raypos += 4 + numsegsperray[sourcecnt]*lSizeofRaypath;
		
		// increment the ray counter
		sourcecnt ++;
	} // for each launched ray

	//mexPrintf("number of source rays = %d, number of segments found = %d, last source object = %d\n", sourcecnt, segsfoundcnt, src_object[sourcecnt-1]);
    
	// calculate Irradiance and field values from rays for each source
	ProcessRays_Parms sParms;
	sParms.Nsegs = segsfoundcnt;
	sParms.bufferpos_segs = bufferpos_rayseg;
	sParms.pdWavelength = pdWavelength;
	sParms.dDetectorSizeX = pdDetectorSize[0];
	sParms.dDetectorSizeY = pdDetectorSize[1];
	sParms.NPixX = (int)(pdNPix[0]);
	sParms.NPixY = (int)(pdNPix[1]);
	sParms.pdIncIrr = (double*)mxCalloc(sParms.NPixX * sParms.NPixY, sizeof(double));
	sParms.pdFieldRe = (double*)mxCalloc(sParms.NPixX * sParms.NPixY, sizeof(double));
	sParms.pdFieldIm = (double*)mxCalloc(sParms.NPixX * sParms.NPixY, sizeof(double));
	sParms.pdSumAmp = (double*)mxCalloc(sParms.NPixX * sParms.NPixY, sizeof(double));
	sParms.px = (double*)mxCalloc(sParms.NPixX, sizeof(double));
	sParms.py = (double*)mxCalloc(sParms.NPixY, sizeof(double));

	ProcessRays(&sParms);

	// send detector arrays outputs to plhs[0]
	OutputMapsToStruct(&sParms, &(plhs[0]));

	// output struct array of found ray segments, if requested
	if (bReturnRays)
		Segment2mxStruct(segsfoundcnt, bufferpos_rayseg, pdWavelength, &(plhs[1]));
	
	////CopyArrayToOutput(sourcecnt, num_dethits, &(plhs[0]));

	mxFree(pSIT);

	mxFree(numsegsperray);
	mxFree(src_object);

	mxFree(pdWavelength);
    
	mxFree(bufferpos_rayseg);

} // mexFunction
