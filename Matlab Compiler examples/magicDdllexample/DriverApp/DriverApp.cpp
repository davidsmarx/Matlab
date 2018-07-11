// DriverApp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "magicsquare.h"

int _tmain(int argc, _TCHAR* argv[])
{

	bool bVal = mclInitializeApplication(NULL,0);
	printf("bVal = %d\n",bVal);
	//if (!libmatrixInitialize()) {
	//	printf("could not initialize the library!\n");
	//	return -1;
	//}
	
	bVal = magicsquareInitialize();


	//bool mlfMagicsquare(int nargout, mxArray** m, mxArray* n);
	mxArray* mm = NULL;
	mxArray* nn;
	nn = mxCreateDoubleMatrix(1,1,mxREAL);
	double nsize = 3.0;
	memcpy(mxGetPr(nn),&nsize,1*sizeof(double));

	bVal = mlfMagicsquare(1, &mm, nn);

	return 0;
}

