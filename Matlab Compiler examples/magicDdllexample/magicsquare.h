/*
 * MATLAB Compiler: 4.6 (R2007a)
 * Date: Thu Aug 09 10:30:09 2007
 * Arguments: "-B" "macro_default" "-l" "magicsquare.m" 
 */

#ifndef __magicsquare_h
#define __magicsquare_h 1

#if defined(__cplusplus) && !defined(mclmcr_h) && defined(__linux__)
#  pragma implementation "mclmcr.h"
#endif
#include "mclmcr.h"
#ifdef __cplusplus
extern "C" {
#endif

#if defined(__SUNPRO_CC)
/* Solaris shared libraries use __global, rather than mapfiles
 * to define the API exported from a shared library. __global is
 * only necessary when building the library -- files including
 * this header file to use the library do not need the __global
 * declaration; hence the EXPORTING_<library> logic.
 */

#ifdef EXPORTING_magicsquare
#define PUBLIC_magicsquare_C_API __global
#else
#define PUBLIC_magicsquare_C_API /* No import statement needed. */
#endif

#define LIB_magicsquare_C_API PUBLIC_magicsquare_C_API

#elif defined(_HPUX_SOURCE)

#ifdef EXPORTING_magicsquare
#define PUBLIC_magicsquare_C_API __declspec(dllexport)
#else
#define PUBLIC_magicsquare_C_API __declspec(dllimport)
#endif

#define LIB_magicsquare_C_API PUBLIC_magicsquare_C_API


#else

#define LIB_magicsquare_C_API

#endif

/* This symbol is defined in shared libraries. Define it here
 * (to nothing) in case this isn't a shared library. 
 */
#ifndef LIB_magicsquare_C_API 
#define LIB_magicsquare_C_API /* No special import/export declaration */
#endif

extern LIB_magicsquare_C_API 
bool MW_CALL_CONV magicsquareInitializeWithHandlers(mclOutputHandlerFcn error_handler,
                                                    mclOutputHandlerFcn print_handler);

extern LIB_magicsquare_C_API 
bool MW_CALL_CONV magicsquareInitialize(void);

extern LIB_magicsquare_C_API 
void MW_CALL_CONV magicsquareTerminate(void);


extern LIB_magicsquare_C_API 
bool MW_CALL_CONV mlxMagicsquare(int nlhs, mxArray *plhs[],
                                 int nrhs, mxArray *prhs[]);


extern LIB_magicsquare_C_API bool MW_CALL_CONV mlfMagicsquare(int nargout
                                                              , mxArray** m
                                                              , mxArray* n);

#ifdef __cplusplus
}
#endif

#endif
