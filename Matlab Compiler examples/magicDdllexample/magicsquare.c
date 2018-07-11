/*
 * MATLAB Compiler: 4.6 (R2007a)
 * Date: Thu Aug 09 10:30:09 2007
 * Arguments: "-B" "macro_default" "-l" "magicsquare.m" 
 */

#include <stdio.h>
#define EXPORTING_magicsquare 1
#include "magicsquare.h"
#ifdef __cplusplus
extern "C" {
#endif

extern mclComponentData __MCC_magicsquare_component_data;

#ifdef __cplusplus
}
#endif


static HMCRINSTANCE _mcr_inst = NULL;


#if defined( _MSC_VER) || defined(__BORLANDC__) || defined(__WATCOMC__) || defined(__LCC__)
#ifdef __LCC__
#undef EXTERN_C
#endif
#include <windows.h>

static char path_to_dll[_MAX_PATH];

BOOL WINAPI DllMain(HINSTANCE hInstance, DWORD dwReason, void *pv)
{
    if (dwReason == DLL_PROCESS_ATTACH)
    {
        char szDllPath[_MAX_PATH];
        char szDir[_MAX_DIR];
        if (GetModuleFileName(hInstance, szDllPath, _MAX_PATH) > 0)
        {
             _splitpath(szDllPath, path_to_dll, szDir, NULL, NULL);
            strcat(path_to_dll, szDir);
        }
	else return FALSE;
    }
    else if (dwReason == DLL_PROCESS_DETACH)
    {
    }
    return TRUE;
}
#endif
#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultPrintHandler(const char *s)
{
    return mclWrite(1 /* stdout */, s, sizeof(char)*strlen(s));
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

#ifdef __cplusplus
extern "C" {
#endif

static int mclDefaultErrorHandler(const char *s)
{
    int written = 0;
    size_t len = 0;
    len = strlen(s);
    written = mclWrite(2 /* stderr */, s, sizeof(char)*len);
    if (len > 0 && s[ len-1 ] != '\n')
        written += mclWrite(2 /* stderr */, "\n", sizeof(char));
    return written;
}

#ifdef __cplusplus
} /* End extern "C" block */
#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_magicsquare_C_API 
#define LIB_magicsquare_C_API /* No special import/export declaration */
#endif

LIB_magicsquare_C_API 
bool MW_CALL_CONV magicsquareInitializeWithHandlers(
    mclOutputHandlerFcn error_handler,
    mclOutputHandlerFcn print_handler
)
{
    if (_mcr_inst != NULL)
        return true;
    if (!mclmcrInitialize())
        return false;
    if (!mclInitializeComponentInstance(&_mcr_inst,
                                        &__MCC_magicsquare_component_data,
                                        true, NoObjectType, LibTarget,
                                        error_handler, print_handler))
        return false;
    return true;
}

LIB_magicsquare_C_API 
bool MW_CALL_CONV magicsquareInitialize(void)
{
    return magicsquareInitializeWithHandlers(mclDefaultErrorHandler,
                                             mclDefaultPrintHandler);
}

LIB_magicsquare_C_API 
void MW_CALL_CONV magicsquareTerminate(void)
{
    if (_mcr_inst != NULL)
        mclTerminateInstance(&_mcr_inst);
}


LIB_magicsquare_C_API 
bool MW_CALL_CONV mlxMagicsquare(int nlhs, mxArray *plhs[],
                                 int nrhs, mxArray *prhs[])
{
    return mclFeval(_mcr_inst, "magicsquare", nlhs, plhs, nrhs, prhs);
}

LIB_magicsquare_C_API 
bool MW_CALL_CONV mlfMagicsquare(int nargout, mxArray** m, mxArray* n)
{
    return mclMlfFeval(_mcr_inst, "magicsquare", nargout, 1, 1, m, n);
}
