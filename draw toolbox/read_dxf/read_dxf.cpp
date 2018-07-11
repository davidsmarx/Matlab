////////////////////////////////////////////////////////////
//
// Name:    readdxf.h
//
// Author:  Steven Michael
//          (smichael@ll.mit.edu)
//
// Date:    3/10/2005
//
// Description:
//
//   The following implements a MATLAB interface to the
//   DXF C++ class that allows for reading of DXF files
//   into MATLAB.
//
////////////////////////////////////////////////////////////

#include <mex.h>

#include <string.h>

#include "dxf.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  if(nrhs<1) {
    mexPrintf("First argument must be a filename.\n");
    return;
  }
  if(!mxIsChar(prhs[0])) {
    mexPrintf("First argument must be a filename.\n");
    return;
  }
  char filename[256];
  mxGetString(prhs[0],filename,256);
  
  DXF *dxf = new DXF;
  if(dxf->read_file(filename)) {
    mexPrintf("Error reading DXF file: \"%s\"\n",
	      filename);
    return;
  }
  if(dxf->load_facets()) {
    mexPrintf("Error loading facets from \"%s\"\n",
	      filename);
    return;
  }
  
  int dim[3];
  dim[2] = dxf->nFacet;
  dim[0] = 3;
  dim[1] = 3;
  
  plhs[0] = mxCreateNumericArray(3,dim,mxSINGLE_CLASS,mxREAL);
  float *data = (float *) mxGetPr(plhs[0]);

  memcpy(data,(const void *)dxf->facets,
	 sizeof(float)*dim[0]*dim[1]*dim[2]);

  delete dxf;

  return;
} // end of mexFunction
