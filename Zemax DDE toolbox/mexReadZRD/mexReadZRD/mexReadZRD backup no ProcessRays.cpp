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

#include "ReadBufferCalcFields.h"

using namespace std;

#define MAXNUMTHREADS 22

class CLogFile
{
	FILE *m_logFile;
	bool m_bLogginOn;

public:
	CLogFile(char *filename)
	{
		m_bLogginOn = false;
		if ( !m_bLogginOn ) return;
		fopen_s(&m_logFile, filename, "at");
	}
	~CLogFile()
	{
		if ( !m_bLogginOn ) return;
		fclose(m_logFile);
	}
	void Append(char *msg)
	{
		if ( !m_bLogginOn ) return;
		fprintf(m_logFile, msg);
		fflush(m_logFile);
	}
	void Append(char *msg, long lNum)
	{
		if ( !m_bLogginOn ) return;
		fprintf(m_logFile, msg, lNum);
		fflush(m_logFile);
	}
};


typedef struct	
{
	size_t N;
	INT32* pList;
} SOURCELIST;

void CopyTranspose(double *pIn, double *pOut, int Nx, int Ny)
{
	for ( int ix = 0; ix < Nx; ix++ )
		for ( int iy = 0; iy < Ny; iy++ )
			pOut[iy + Ny*ix] = pIn[ix + Nx*iy];
} // CopyTranspose

void OutputMapsToStruct(ReadBufferCalcFields_Parms *pSMaps, SOURCELIST *psSources,  mxArray **pmxMapStruct)
{	
	const char *field_names[] = { "IncIrr", "Field", "SumAmp", "x", "y", "SourceObj" };
	const long lNumFields = (sizeof(field_names)/sizeof(*field_names));

	// create struct 1x1 array
	mwSize dims[2] = {psSources->N, 1};
	*pmxMapStruct = mxCreateStructArray(2, dims, lNumFields, field_names);
	/* Since we just
       created the structure and the field number indices are zero
       based, , etc. */
	
	for ( int isou = 0; isou < psSources->N; isou++ )
	{
		mxArray *IncIrr = mxCreateDoubleMatrix(pSMaps[isou].NPixX, pSMaps[isou].NPixY, mxREAL);
		double *pointer = mxGetPr(IncIrr);
	   // memcpy(pointer, pSMaps->pdIncIrr, pSMaps->NPixX * pSMaps->NPixY * sizeof(double));
		CopyTranspose(pSMaps[isou].pdIncIrr, pointer, pSMaps[isou].NPixX, pSMaps[isou].NPixY);
		mxSetFieldByNumber(*pmxMapStruct, isou, 0, IncIrr);

		mxArray *CohField = mxCreateDoubleMatrix(pSMaps[isou].NPixX, pSMaps[isou].NPixY, mxCOMPLEX);
		double *pRe = mxGetPr(CohField);
		double *pIm = mxGetPi(CohField);
		//memcpy(pRe, pSMaps->pdFieldRe,  pSMaps->NPixX * pSMaps->NPixY * sizeof(double));
		CopyTranspose(pSMaps[isou].pdFieldRe, pRe, pSMaps[isou].NPixX, pSMaps[isou].NPixY);
		//memcpy(pIm, pSMaps->pdFieldIm,  pSMaps->NPixX * pSMaps->NPixY * sizeof(double));
		CopyTranspose(pSMaps[isou].pdFieldIm, pIm, pSMaps[isou].NPixX, pSMaps[isou].NPixY);
		mxSetFieldByNumber(*pmxMapStruct, isou, 1, CohField);

		mxArray *SumAmp = mxCreateDoubleMatrix(pSMaps[isou].NPixX, pSMaps[isou].NPixY, mxREAL);
		double *pSA = mxGetPr(SumAmp);
		CopyTranspose(pSMaps[isou].pdSumAmp, pSA, pSMaps[isou].NPixX, pSMaps[isou].NPixY);
		mxSetFieldByNumber(*pmxMapStruct, isou, 2, SumAmp);

		mxArray *x = mxCreateDoubleMatrix(pSMaps[isou].NPixX, 1, mxREAL );
		double *px = mxGetPr(x);
		memcpy(px, pSMaps[isou].px, pSMaps[isou].NPixX * sizeof(double));
		mxSetFieldByNumber(*pmxMapStruct, isou, 3, x);

		mxArray *y = mxCreateDoubleMatrix(pSMaps[isou].NPixY, 1, mxREAL );
		double *py = mxGetPr(y);
		memcpy(py, pSMaps[isou].py, pSMaps[isou].NPixY * sizeof(double));
		mxSetFieldByNumber(*pmxMapStruct, isou, 4, y);

		mxArray *SourceObj = mxCreateDoubleMatrix(1, 1, mxREAL );
		double *pS = mxGetPr(SourceObj);
		*pS = psSources->pList[isou];
		mxSetFieldByNumber(*pmxMapStruct, isou, 5, SourceObj);

	} // for each source

} // OutputMapsToStruct

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    
	CLogFile LogFile("C:\\Users\\davidmarx\\Documents\\multi incoh sources 2014-05-02\\test all sources to zrd\\logfile.txt");

	LogFile.Append("start mexFunction\n");
    
    /* Check for proper number of input and  output arguments */    
    if (nrhs != 2) {
        mexErrMsgIdAndTxt( "MATLAB:mxcreatestructarray:maxrhs",
                "usage: mexReadZRD(filename, sOptions)");
    } 


	if ( nlhs != 1 ) {
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
	SOURCELIST sSourceList;

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


	// calculate Irradiance and field values from rays for each source
	ReadBufferCalcFields_Parms *sParms;
	sParms = (ReadBufferCalcFields_Parms*)mxCalloc(sSourceList.N, sizeof(ReadBufferCalcFields_Parms));

	HANDLE *pHandle = new HANDLE[MAXNUMTHREADS];
	int cntthrd = 0;
	for ( int isou = 0; isou < sSourceList.N; isou++ )
	{
		int Nx = (int)(pdNPix[0]);
		int Ny = (int)(pdNPix[1]);
		sParms[isou].iSourceObj = sSourceList.pList[isou];
		sParms[isou].dDetectorSizeX = pdDetectorSize[0];
		sParms[isou].dDetectorSizeY = pdDetectorSize[1];
		sParms[isou].NPixX = Nx;
		sParms[isou].NPixY = Ny;
		sParms[isou].buffersize = buffersize;
		sParms[isou].buffer = buffer;
		sParms[isou].lFirstSegment = lFirstSegment;
		sParms[isou].lDetectorObject = lDetectorObject;

		sParms[isou].pdIncIrr = (double*)mxCalloc(Nx * Ny, sizeof(double));
		sParms[isou].pdFieldRe = (double*)mxCalloc(Nx * Ny, sizeof(double));
		sParms[isou].pdFieldIm = (double*)mxCalloc(Nx * Ny, sizeof(double));
		sParms[isou].pdSumAmp = (double*)mxCalloc(Nx * Ny, sizeof(double));
		sParms[isou].px = (double*)mxCalloc(Nx, sizeof(double));
		sParms[isou].py = (double*)mxCalloc(Ny, sizeof(double));
		
		sParms[isou].pNumSegsFound = (int*)mxCalloc(1, sizeof(int));
	
		//ReadBufferCalcFields(&(sParms[isou]));
		pHandle[cntthrd] = CreateThread( 
			NULL,						// no security attributes 
			0,							// use default stack size  
			ReadBufferCalcFields,// thread function to wait for scan completion
			&(sParms[isou]),   	 // argument to thread function 
			0,							// use default creation flags 
			NULL);						// returns the thread identifier 

        if (pHandle[cntthrd] == NULL)
		{
			mexPrintf("CreateThread failed for isou = %d, Source Obj = %d\n",isou, sSourceList.pList[isou]);
			mexErrMsgIdAndTxt( "mexReadZRD:CreateThread", "CreateThread Failed");
		}

		cntthrd++;
		
		LogFile.Append("created thread # %d", cntthrd);
		LogFile.Append(" for source number %d\n", sSourceList.pList[isou]);

		if ( cntthrd == MAXNUMTHREADS )
		{
			// wait for thread to terminate
			WaitForMultipleObjects(cntthrd, pHandle, true, INFINITE);

			// recycle the handles
			for ( int ithr = 0; ithr < MAXNUMTHREADS; ithr++ )
				CloseHandle(pHandle[ithr]);

			LogFile.Append("wait for threads done, setting cntthrd = 0\n");

			cntthrd = 0;
		}

	} // for each source object

	// wait for all the remaining threads to finish
	WaitForMultipleObjects(cntthrd, pHandle, true, INFINITE);

	for ( int ithr = 0; ithr < cntthrd; ithr++ )
		CloseHandle(pHandle[ithr]);

	LogFile.Append("done waiting for last %d threads\n",cntthrd);

	// send detector arrays outputs to plhs[0]
	OutputMapsToStruct(sParms, &sSourceList, &(plhs[0]));
	
	LogFile.Append("Done OutputMapsToStruct...");

	delete [] pHandle;

	LogFile.Append("Done everything.\n");

} // mexFunction
