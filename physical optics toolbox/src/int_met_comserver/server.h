/* this ALWAYS GENERATED file contains the definitions for the interfaces */


/* File created by MIDL compiler version 5.01.0164 */
/* at Tue Aug 09 12:28:40 2005
 */
/* Compiler settings for C:\My_Documents\internal_metrology\diffraction\int_met_comserver\server.idl:
    Os (OptLev=s), W1, Zp8, env=Win32, ms_ext, c_ext
    error checks: allocation ref bounds_check enum stub_data 
*/
//@@MIDL_FILE_HEADING(  )


/* verify that the <rpcndr.h> version is high enough to compile this file*/
#ifndef __REQUIRED_RPCNDR_H_VERSION__
#define __REQUIRED_RPCNDR_H_VERSION__ 440
#endif

#include "rpc.h"
#include "rpcndr.h"

#ifndef __RPCNDR_H_VERSION__
#error this stub requires an updated version of <rpcndr.h>
#endif // __RPCNDR_H_VERSION__

#ifndef COM_NO_WINDOWS_H
#include "windows.h"
#include "ole2.h"
#endif /*COM_NO_WINDOWS_H*/

#ifndef __server_h__
#define __server_h__

#ifdef __cplusplus
extern "C"{
#endif 

/* Forward Declarations */ 

#ifndef __IDiffraction_FWD_DEFINED__
#define __IDiffraction_FWD_DEFINED__
typedef interface IDiffraction IDiffraction;
#endif 	/* __IDiffraction_FWD_DEFINED__ */


#ifndef __Diffraction_FWD_DEFINED__
#define __Diffraction_FWD_DEFINED__

#ifdef __cplusplus
typedef class Diffraction Diffraction;
#else
typedef struct Diffraction Diffraction;
#endif /* __cplusplus */

#endif 	/* __Diffraction_FWD_DEFINED__ */


/* header files for imported files */
#include "oaidl.h"
#include "ocidl.h"

void __RPC_FAR * __RPC_USER MIDL_user_allocate(size_t);
void __RPC_USER MIDL_user_free( void __RPC_FAR * ); 

#ifndef __IDiffraction_INTERFACE_DEFINED__
#define __IDiffraction_INTERFACE_DEFINED__

/* interface IDiffraction */
/* [dual][oleautomation][uuid][object] */ 


EXTERN_C const IID IID_IDiffraction;

#if defined(__cplusplus) && !defined(CINTERFACE)
    
    MIDL_INTERFACE("44182096-CE96-4DC8-8D63-EC2FC6FA8D85")
    IDiffraction : public IDispatch
    {
    public:
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE x_vector_get( 
            /* [in] */ long BeamID,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *x_vector) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE y_vector_get( 
            /* [in] */ long BeamID,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *y_vector) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE wavefront_put( 
            /* [in] */ SAFEARRAY __RPC_FAR * amp_r,
            /* [in] */ SAFEARRAY __RPC_FAR * amp_i,
            /* [in] */ double dx,
            /* [in] */ double dy,
            /* [in] */ double wavelength,
            /* [in] */ double curv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE wavefront_get( 
            /* [in] */ long BeamID,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *amp_r,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *amp_i) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE CopyWavefront( 
            /* [in] */ long BeamIDfrom,
            /* [in] */ long BeamIDto,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE CopyWavefrontAmp( 
            /* [in] */ long BeamIDAmp,
            /* [in] */ long BeamIDPhase,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE SumWavefront( 
            /* [in] */ long BeamID,
            /* [in] */ long BeamIDadd,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE GetWavefrontParms( 
            /* [in] */ long BeamID,
            /* [out] */ long __RPC_FAR *Nx,
            /* [out] */ long __RPC_FAR *Ny,
            /* [out] */ double __RPC_FAR *dx,
            /* [out] */ double __RPC_FAR *dy,
            /* [out] */ double __RPC_FAR *curv,
            /* [out] */ double __RPC_FAR *wavelength) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE WriteWavefrontUNF( 
            /* [in] */ BSTR StrFilename,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ReadWavefrontUNF( 
            /* [in] */ BSTR StrFilename,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE WavefrontPower( 
            /* [out] */ double __RPC_FAR *power,
            /* [in] */ long BeamID) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE CreateGausSource( 
            /* [in] */ double BeamWaistDiam,
            /* [in] */ long Nx,
            /* [in] */ long Ny,
            /* [in] */ double dx,
            /* [in] */ double dy,
            /* [in] */ double wavelength,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE CreateTophatSource( 
            /* [in] */ double BeamDiam,
            /* [in] */ long Nx,
            /* [in] */ long Ny,
            /* [in] */ double dx,
            /* [in] */ double dy,
            /* [in] */ double wavelength,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ClearWavefront( 
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMask( 
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMaskMisalign( 
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMaskGeneral( 
            /* [in] */ double len_x,
            /* [in] */ double len_y,
            /* [in] */ double off_x,
            /* [in] */ double off_y,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMaskPoly( 
            /* [in] */ SAFEARRAY __RPC_FAR * vertices,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMaskX( 
            /* [in] */ double radius,
            /* [in] */ double length,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMaskXRounded( 
            /* [in] */ double radius,
            /* [in] */ double length,
            /* [in] */ double corner_rad,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMaskRotate( 
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ double angle,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ApplyMaskTilt( 
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ double rotangle,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double tiltangle,
            /* [in] */ double tiltorientation,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE Propagate( 
            /* [in] */ double Distance,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE PropagateExt( 
            /* [in] */ double Distance,
            /* [in] */ double dxout,
            /* [in] */ double dyout,
            /* [in] */ long applyCurv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE PropagateDefDx( 
            /* [in] */ double Distance,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE CornerCube( 
            /* [in] */ double ccsize,
            /* [in] */ BSTR shape,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double spin,
            /* [in] */ SAFEARRAY __RPC_FAR * gapwidth,
            /* [in] */ SAFEARRAY __RPC_FAR * rotmatrix,
            /* [in] */ SAFEARRAY __RPC_FAR * dihedral,
            /* [in] */ SAFEARRAY __RPC_FAR * edgelength,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE focuslens( 
            /* [in] */ double focuslens_f,
            /* [in] */ double focuslens_D,
            /* [in] */ double dxout,
            /* [in] */ double dyout,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ThinLens( 
            /* [in] */ double focallength,
            /* [in] */ double diameter,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ReflectAsphere( 
            /* [in] */ double rCurv,
            /* [in] */ double conicConst,
            /* [in] */ double diam,
            /* [in] */ double xDecenter,
            /* [in] */ double yDecenter,
            /* [in] */ double incidenceAngle,
            /* [in] */ double azimuth,
            /* [in] */ long applyCurv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE PistonTiltFocus( 
            /* [in] */ double piston,
            /* [in] */ double xtilt,
            /* [in] */ double ytilt,
            /* [in] */ double focus,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE PropagateToFocalPlane( 
            /* [in] */ double lens_f,
            /* [in] */ double lens_D,
            /* [in] */ double dxout,
            /* [in] */ double dyout,
            /* [in] */ long applycurv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ClipCirc( 
            /* [in] */ double Diameter,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double obscDiam,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ClipRect( 
            /* [in] */ double len_x,
            /* [in] */ double len_y,
            /* [in] */ double off_x,
            /* [in] */ double off_y,
            /* [in] */ long obscFlag,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ClipPolygon( 
            /* [in] */ SAFEARRAY __RPC_FAR * vertices,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double angle,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ClipBWindow( 
            /* [in] */ double winT,
            /* [in] */ double winA,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double rotAngle,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE AddThreadCommand( 
            /* [in] */ long ThreadID,
            /* [in] */ long CmdID,
            /* [in] */ SAFEARRAY __RPC_FAR * rP,
            /* [in] */ SAFEARRAY __RPC_FAR * iP,
            /* [in] */ long pwAID,
            /* [in] */ long pwBID,
            /* [in] */ BSTR strP,
            /* [retval][out] */ long __RPC_FAR *flagStatus) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ExecuteThread( 
            /* [in] */ long ThreadID,
            /* [retval][out] */ long __RPC_FAR *flagStatus) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ThreadStatus( 
            /* [in] */ long ThreadID,
            /* [retval][out] */ long __RPC_FAR *flagStatus) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE WaitForThread( 
            /* [in] */ long ThreadID,
            /* [in] */ long TimeOutMS,
            /* [retval][out] */ long __RPC_FAR *WaitStatus) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE WaitForMultipleThreads( 
            /* [in] */ SAFEARRAY __RPC_FAR * pThreadID,
            /* [in] */ long TimeOutMS,
            /* [in] */ long flagWaitAll,
            /* [retval][out] */ long __RPC_FAR *WaitStatus) = 0;
        
        virtual /* [id] */ HRESULT STDMETHODCALLTYPE ClearThread( 
            /* [in] */ long ThreadID,
            /* [retval][out] */ long __RPC_FAR *flagStatus) = 0;
        
        virtual /* [propget][id] */ HRESULT STDMETHODCALLTYPE get_ThreadOPD( 
            /* [in] */ long ThreadID,
            /* [in] */ long OPDID,
            /* [retval][out] */ double __RPC_FAR *VALUE) = 0;
        
        virtual /* [propget][id] */ HRESULT STDMETHODCALLTYPE get_ThreadPOW( 
            /* [in] */ long ThreadID,
            /* [in] */ long POWID,
            /* [retval][out] */ double __RPC_FAR *VALUE) = 0;
        
    };
    
#else 	/* C style interface */

    typedef struct IDiffractionVtbl
    {
        BEGIN_INTERFACE
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *QueryInterface )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [iid_is][out] */ void __RPC_FAR *__RPC_FAR *ppvObject);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *AddRef )( 
            IDiffraction __RPC_FAR * This);
        
        ULONG ( STDMETHODCALLTYPE __RPC_FAR *Release )( 
            IDiffraction __RPC_FAR * This);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTypeInfoCount )( 
            IDiffraction __RPC_FAR * This,
            /* [out] */ UINT __RPC_FAR *pctinfo);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetTypeInfo )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ UINT iTInfo,
            /* [in] */ LCID lcid,
            /* [out] */ ITypeInfo __RPC_FAR *__RPC_FAR *ppTInfo);
        
        HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetIDsOfNames )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ REFIID riid,
            /* [size_is][in] */ LPOLESTR __RPC_FAR *rgszNames,
            /* [in] */ UINT cNames,
            /* [in] */ LCID lcid,
            /* [size_is][out] */ DISPID __RPC_FAR *rgDispId);
        
        /* [local] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Invoke )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ DISPID dispIdMember,
            /* [in] */ REFIID riid,
            /* [in] */ LCID lcid,
            /* [in] */ WORD wFlags,
            /* [out][in] */ DISPPARAMS __RPC_FAR *pDispParams,
            /* [out] */ VARIANT __RPC_FAR *pVarResult,
            /* [out] */ EXCEPINFO __RPC_FAR *pExcepInfo,
            /* [out] */ UINT __RPC_FAR *puArgErr);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *x_vector_get )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamID,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *x_vector);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *y_vector_get )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamID,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *y_vector);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *wavefront_put )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ SAFEARRAY __RPC_FAR * amp_r,
            /* [in] */ SAFEARRAY __RPC_FAR * amp_i,
            /* [in] */ double dx,
            /* [in] */ double dy,
            /* [in] */ double wavelength,
            /* [in] */ double curv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *wavefront_get )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamID,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *amp_r,
            /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *amp_i);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *CopyWavefront )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamIDfrom,
            /* [in] */ long BeamIDto,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *CopyWavefrontAmp )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamIDAmp,
            /* [in] */ long BeamIDPhase,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *SumWavefront )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamID,
            /* [in] */ long BeamIDadd,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *GetWavefrontParms )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamID,
            /* [out] */ long __RPC_FAR *Nx,
            /* [out] */ long __RPC_FAR *Ny,
            /* [out] */ double __RPC_FAR *dx,
            /* [out] */ double __RPC_FAR *dy,
            /* [out] */ double __RPC_FAR *curv,
            /* [out] */ double __RPC_FAR *wavelength);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *WriteWavefrontUNF )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ BSTR StrFilename,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ReadWavefrontUNF )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ BSTR StrFilename,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *WavefrontPower )( 
            IDiffraction __RPC_FAR * This,
            /* [out] */ double __RPC_FAR *power,
            /* [in] */ long BeamID);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *CreateGausSource )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double BeamWaistDiam,
            /* [in] */ long Nx,
            /* [in] */ long Ny,
            /* [in] */ double dx,
            /* [in] */ double dy,
            /* [in] */ double wavelength,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *CreateTophatSource )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double BeamDiam,
            /* [in] */ long Nx,
            /* [in] */ long Ny,
            /* [in] */ double dx,
            /* [in] */ double dy,
            /* [in] */ double wavelength,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ClearWavefront )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMask )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMaskMisalign )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMaskGeneral )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double len_x,
            /* [in] */ double len_y,
            /* [in] */ double off_x,
            /* [in] */ double off_y,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMaskPoly )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ SAFEARRAY __RPC_FAR * vertices,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMaskX )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double radius,
            /* [in] */ double length,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMaskXRounded )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double radius,
            /* [in] */ double length,
            /* [in] */ double corner_rad,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMaskRotate )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ double angle,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ApplyMaskTilt )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double length_r,
            /* [in] */ double length_t,
            /* [in] */ double offset,
            /* [in] */ double rotangle,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double tiltangle,
            /* [in] */ double tiltorientation,
            /* [in] */ BSTR direction,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *Propagate )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double Distance,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PropagateExt )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double Distance,
            /* [in] */ double dxout,
            /* [in] */ double dyout,
            /* [in] */ long applyCurv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PropagateDefDx )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double Distance,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *CornerCube )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double ccsize,
            /* [in] */ BSTR shape,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double spin,
            /* [in] */ SAFEARRAY __RPC_FAR * gapwidth,
            /* [in] */ SAFEARRAY __RPC_FAR * rotmatrix,
            /* [in] */ SAFEARRAY __RPC_FAR * dihedral,
            /* [in] */ SAFEARRAY __RPC_FAR * edgelength,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *focuslens )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double focuslens_f,
            /* [in] */ double focuslens_D,
            /* [in] */ double dxout,
            /* [in] */ double dyout,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ThinLens )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double focallength,
            /* [in] */ double diameter,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ReflectAsphere )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double rCurv,
            /* [in] */ double conicConst,
            /* [in] */ double diam,
            /* [in] */ double xDecenter,
            /* [in] */ double yDecenter,
            /* [in] */ double incidenceAngle,
            /* [in] */ double azimuth,
            /* [in] */ long applyCurv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PistonTiltFocus )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double piston,
            /* [in] */ double xtilt,
            /* [in] */ double ytilt,
            /* [in] */ double focus,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *PropagateToFocalPlane )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double lens_f,
            /* [in] */ double lens_D,
            /* [in] */ double dxout,
            /* [in] */ double dyout,
            /* [in] */ long applycurv,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ClipCirc )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double Diameter,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double obscDiam,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ClipRect )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double len_x,
            /* [in] */ double len_y,
            /* [in] */ double off_x,
            /* [in] */ double off_y,
            /* [in] */ long obscFlag,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ClipPolygon )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ SAFEARRAY __RPC_FAR * vertices,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double angle,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ClipBWindow )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ double winT,
            /* [in] */ double winA,
            /* [in] */ double xc,
            /* [in] */ double yc,
            /* [in] */ double rotAngle,
            /* [in] */ long BeamID,
            /* [retval][out] */ long __RPC_FAR *status);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *AddThreadCommand )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long ThreadID,
            /* [in] */ long CmdID,
            /* [in] */ SAFEARRAY __RPC_FAR * rP,
            /* [in] */ SAFEARRAY __RPC_FAR * iP,
            /* [in] */ long pwAID,
            /* [in] */ long pwBID,
            /* [in] */ BSTR strP,
            /* [retval][out] */ long __RPC_FAR *flagStatus);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ExecuteThread )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long ThreadID,
            /* [retval][out] */ long __RPC_FAR *flagStatus);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ThreadStatus )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long ThreadID,
            /* [retval][out] */ long __RPC_FAR *flagStatus);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *WaitForThread )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long ThreadID,
            /* [in] */ long TimeOutMS,
            /* [retval][out] */ long __RPC_FAR *WaitStatus);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *WaitForMultipleThreads )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ SAFEARRAY __RPC_FAR * pThreadID,
            /* [in] */ long TimeOutMS,
            /* [in] */ long flagWaitAll,
            /* [retval][out] */ long __RPC_FAR *WaitStatus);
        
        /* [id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *ClearThread )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long ThreadID,
            /* [retval][out] */ long __RPC_FAR *flagStatus);
        
        /* [propget][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *get_ThreadOPD )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long ThreadID,
            /* [in] */ long OPDID,
            /* [retval][out] */ double __RPC_FAR *VALUE);
        
        /* [propget][id] */ HRESULT ( STDMETHODCALLTYPE __RPC_FAR *get_ThreadPOW )( 
            IDiffraction __RPC_FAR * This,
            /* [in] */ long ThreadID,
            /* [in] */ long POWID,
            /* [retval][out] */ double __RPC_FAR *VALUE);
        
        END_INTERFACE
    } IDiffractionVtbl;

    interface IDiffraction
    {
        CONST_VTBL struct IDiffractionVtbl __RPC_FAR *lpVtbl;
    };

    

#ifdef COBJMACROS


#define IDiffraction_QueryInterface(This,riid,ppvObject)	\
    (This)->lpVtbl -> QueryInterface(This,riid,ppvObject)

#define IDiffraction_AddRef(This)	\
    (This)->lpVtbl -> AddRef(This)

#define IDiffraction_Release(This)	\
    (This)->lpVtbl -> Release(This)


#define IDiffraction_GetTypeInfoCount(This,pctinfo)	\
    (This)->lpVtbl -> GetTypeInfoCount(This,pctinfo)

#define IDiffraction_GetTypeInfo(This,iTInfo,lcid,ppTInfo)	\
    (This)->lpVtbl -> GetTypeInfo(This,iTInfo,lcid,ppTInfo)

#define IDiffraction_GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)	\
    (This)->lpVtbl -> GetIDsOfNames(This,riid,rgszNames,cNames,lcid,rgDispId)

#define IDiffraction_Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)	\
    (This)->lpVtbl -> Invoke(This,dispIdMember,riid,lcid,wFlags,pDispParams,pVarResult,pExcepInfo,puArgErr)


#define IDiffraction_x_vector_get(This,BeamID,x_vector)	\
    (This)->lpVtbl -> x_vector_get(This,BeamID,x_vector)

#define IDiffraction_y_vector_get(This,BeamID,y_vector)	\
    (This)->lpVtbl -> y_vector_get(This,BeamID,y_vector)

#define IDiffraction_wavefront_put(This,amp_r,amp_i,dx,dy,wavelength,curv,BeamID,status)	\
    (This)->lpVtbl -> wavefront_put(This,amp_r,amp_i,dx,dy,wavelength,curv,BeamID,status)

#define IDiffraction_wavefront_get(This,BeamID,amp_r,amp_i)	\
    (This)->lpVtbl -> wavefront_get(This,BeamID,amp_r,amp_i)

#define IDiffraction_CopyWavefront(This,BeamIDfrom,BeamIDto,status)	\
    (This)->lpVtbl -> CopyWavefront(This,BeamIDfrom,BeamIDto,status)

#define IDiffraction_CopyWavefrontAmp(This,BeamIDAmp,BeamIDPhase,status)	\
    (This)->lpVtbl -> CopyWavefrontAmp(This,BeamIDAmp,BeamIDPhase,status)

#define IDiffraction_SumWavefront(This,BeamID,BeamIDadd,status)	\
    (This)->lpVtbl -> SumWavefront(This,BeamID,BeamIDadd,status)

#define IDiffraction_GetWavefrontParms(This,BeamID,Nx,Ny,dx,dy,curv,wavelength)	\
    (This)->lpVtbl -> GetWavefrontParms(This,BeamID,Nx,Ny,dx,dy,curv,wavelength)

#define IDiffraction_WriteWavefrontUNF(This,StrFilename,BeamID,status)	\
    (This)->lpVtbl -> WriteWavefrontUNF(This,StrFilename,BeamID,status)

#define IDiffraction_ReadWavefrontUNF(This,StrFilename,BeamID,status)	\
    (This)->lpVtbl -> ReadWavefrontUNF(This,StrFilename,BeamID,status)

#define IDiffraction_WavefrontPower(This,power,BeamID)	\
    (This)->lpVtbl -> WavefrontPower(This,power,BeamID)

#define IDiffraction_CreateGausSource(This,BeamWaistDiam,Nx,Ny,dx,dy,wavelength,BeamID,status)	\
    (This)->lpVtbl -> CreateGausSource(This,BeamWaistDiam,Nx,Ny,dx,dy,wavelength,BeamID,status)

#define IDiffraction_CreateTophatSource(This,BeamDiam,Nx,Ny,dx,dy,wavelength,BeamID,status)	\
    (This)->lpVtbl -> CreateTophatSource(This,BeamDiam,Nx,Ny,dx,dy,wavelength,BeamID,status)

#define IDiffraction_ClearWavefront(This,BeamID,status)	\
    (This)->lpVtbl -> ClearWavefront(This,BeamID,status)

#define IDiffraction_ApplyMask(This,length_r,length_t,offset,direction,BeamID,status)	\
    (This)->lpVtbl -> ApplyMask(This,length_r,length_t,offset,direction,BeamID,status)

#define IDiffraction_ApplyMaskMisalign(This,length_r,length_t,offset,xc,yc,direction,BeamID,status)	\
    (This)->lpVtbl -> ApplyMaskMisalign(This,length_r,length_t,offset,xc,yc,direction,BeamID,status)

#define IDiffraction_ApplyMaskGeneral(This,len_x,len_y,off_x,off_y,BeamID,status)	\
    (This)->lpVtbl -> ApplyMaskGeneral(This,len_x,len_y,off_x,off_y,BeamID,status)

#define IDiffraction_ApplyMaskPoly(This,vertices,xc,yc,BeamID,status)	\
    (This)->lpVtbl -> ApplyMaskPoly(This,vertices,xc,yc,BeamID,status)

#define IDiffraction_ApplyMaskX(This,radius,length,BeamID,status)	\
    (This)->lpVtbl -> ApplyMaskX(This,radius,length,BeamID,status)

#define IDiffraction_ApplyMaskXRounded(This,radius,length,corner_rad,BeamID,status)	\
    (This)->lpVtbl -> ApplyMaskXRounded(This,radius,length,corner_rad,BeamID,status)

#define IDiffraction_ApplyMaskRotate(This,length_r,length_t,offset,angle,xc,yc,direction,BeamID,status)	\
    (This)->lpVtbl -> ApplyMaskRotate(This,length_r,length_t,offset,angle,xc,yc,direction,BeamID,status)

#define IDiffraction_ApplyMaskTilt(This,length_r,length_t,offset,rotangle,xc,yc,tiltangle,tiltorientation,direction,BeamID,status)	\
    (This)->lpVtbl -> ApplyMaskTilt(This,length_r,length_t,offset,rotangle,xc,yc,tiltangle,tiltorientation,direction,BeamID,status)

#define IDiffraction_Propagate(This,Distance,BeamID,status)	\
    (This)->lpVtbl -> Propagate(This,Distance,BeamID,status)

#define IDiffraction_PropagateExt(This,Distance,dxout,dyout,applyCurv,BeamID,status)	\
    (This)->lpVtbl -> PropagateExt(This,Distance,dxout,dyout,applyCurv,BeamID,status)

#define IDiffraction_PropagateDefDx(This,Distance,BeamID,status)	\
    (This)->lpVtbl -> PropagateDefDx(This,Distance,BeamID,status)

#define IDiffraction_CornerCube(This,ccsize,shape,xc,yc,spin,gapwidth,rotmatrix,dihedral,edgelength,BeamID,status)	\
    (This)->lpVtbl -> CornerCube(This,ccsize,shape,xc,yc,spin,gapwidth,rotmatrix,dihedral,edgelength,BeamID,status)

#define IDiffraction_focuslens(This,focuslens_f,focuslens_D,dxout,dyout,BeamID,status)	\
    (This)->lpVtbl -> focuslens(This,focuslens_f,focuslens_D,dxout,dyout,BeamID,status)

#define IDiffraction_ThinLens(This,focallength,diameter,xc,yc,BeamID,status)	\
    (This)->lpVtbl -> ThinLens(This,focallength,diameter,xc,yc,BeamID,status)

#define IDiffraction_ReflectAsphere(This,rCurv,conicConst,diam,xDecenter,yDecenter,incidenceAngle,azimuth,applyCurv,BeamID,status)	\
    (This)->lpVtbl -> ReflectAsphere(This,rCurv,conicConst,diam,xDecenter,yDecenter,incidenceAngle,azimuth,applyCurv,BeamID,status)

#define IDiffraction_PistonTiltFocus(This,piston,xtilt,ytilt,focus,xc,yc,BeamID,status)	\
    (This)->lpVtbl -> PistonTiltFocus(This,piston,xtilt,ytilt,focus,xc,yc,BeamID,status)

#define IDiffraction_PropagateToFocalPlane(This,lens_f,lens_D,dxout,dyout,applycurv,BeamID,status)	\
    (This)->lpVtbl -> PropagateToFocalPlane(This,lens_f,lens_D,dxout,dyout,applycurv,BeamID,status)

#define IDiffraction_ClipCirc(This,Diameter,xc,yc,obscDiam,BeamID,status)	\
    (This)->lpVtbl -> ClipCirc(This,Diameter,xc,yc,obscDiam,BeamID,status)

#define IDiffraction_ClipRect(This,len_x,len_y,off_x,off_y,obscFlag,BeamID,status)	\
    (This)->lpVtbl -> ClipRect(This,len_x,len_y,off_x,off_y,obscFlag,BeamID,status)

#define IDiffraction_ClipPolygon(This,vertices,xc,yc,angle,BeamID,status)	\
    (This)->lpVtbl -> ClipPolygon(This,vertices,xc,yc,angle,BeamID,status)

#define IDiffraction_ClipBWindow(This,winT,winA,xc,yc,rotAngle,BeamID,status)	\
    (This)->lpVtbl -> ClipBWindow(This,winT,winA,xc,yc,rotAngle,BeamID,status)

#define IDiffraction_AddThreadCommand(This,ThreadID,CmdID,rP,iP,pwAID,pwBID,strP,flagStatus)	\
    (This)->lpVtbl -> AddThreadCommand(This,ThreadID,CmdID,rP,iP,pwAID,pwBID,strP,flagStatus)

#define IDiffraction_ExecuteThread(This,ThreadID,flagStatus)	\
    (This)->lpVtbl -> ExecuteThread(This,ThreadID,flagStatus)

#define IDiffraction_ThreadStatus(This,ThreadID,flagStatus)	\
    (This)->lpVtbl -> ThreadStatus(This,ThreadID,flagStatus)

#define IDiffraction_WaitForThread(This,ThreadID,TimeOutMS,WaitStatus)	\
    (This)->lpVtbl -> WaitForThread(This,ThreadID,TimeOutMS,WaitStatus)

#define IDiffraction_WaitForMultipleThreads(This,pThreadID,TimeOutMS,flagWaitAll,WaitStatus)	\
    (This)->lpVtbl -> WaitForMultipleThreads(This,pThreadID,TimeOutMS,flagWaitAll,WaitStatus)

#define IDiffraction_ClearThread(This,ThreadID,flagStatus)	\
    (This)->lpVtbl -> ClearThread(This,ThreadID,flagStatus)

#define IDiffraction_get_ThreadOPD(This,ThreadID,OPDID,VALUE)	\
    (This)->lpVtbl -> get_ThreadOPD(This,ThreadID,OPDID,VALUE)

#define IDiffraction_get_ThreadPOW(This,ThreadID,POWID,VALUE)	\
    (This)->lpVtbl -> get_ThreadPOW(This,ThreadID,POWID,VALUE)

#endif /* COBJMACROS */


#endif 	/* C style interface */



/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_x_vector_get_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamID,
    /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *x_vector);


void __RPC_STUB IDiffraction_x_vector_get_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_y_vector_get_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamID,
    /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *y_vector);


void __RPC_STUB IDiffraction_y_vector_get_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_wavefront_put_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ SAFEARRAY __RPC_FAR * amp_r,
    /* [in] */ SAFEARRAY __RPC_FAR * amp_i,
    /* [in] */ double dx,
    /* [in] */ double dy,
    /* [in] */ double wavelength,
    /* [in] */ double curv,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_wavefront_put_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_wavefront_get_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamID,
    /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *amp_r,
    /* [out][in] */ SAFEARRAY __RPC_FAR * __RPC_FAR *amp_i);


void __RPC_STUB IDiffraction_wavefront_get_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_CopyWavefront_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamIDfrom,
    /* [in] */ long BeamIDto,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_CopyWavefront_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_CopyWavefrontAmp_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamIDAmp,
    /* [in] */ long BeamIDPhase,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_CopyWavefrontAmp_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_SumWavefront_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamID,
    /* [in] */ long BeamIDadd,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_SumWavefront_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_GetWavefrontParms_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamID,
    /* [out] */ long __RPC_FAR *Nx,
    /* [out] */ long __RPC_FAR *Ny,
    /* [out] */ double __RPC_FAR *dx,
    /* [out] */ double __RPC_FAR *dy,
    /* [out] */ double __RPC_FAR *curv,
    /* [out] */ double __RPC_FAR *wavelength);


void __RPC_STUB IDiffraction_GetWavefrontParms_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_WriteWavefrontUNF_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ BSTR StrFilename,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_WriteWavefrontUNF_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ReadWavefrontUNF_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ BSTR StrFilename,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ReadWavefrontUNF_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_WavefrontPower_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [out] */ double __RPC_FAR *power,
    /* [in] */ long BeamID);


void __RPC_STUB IDiffraction_WavefrontPower_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_CreateGausSource_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double BeamWaistDiam,
    /* [in] */ long Nx,
    /* [in] */ long Ny,
    /* [in] */ double dx,
    /* [in] */ double dy,
    /* [in] */ double wavelength,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_CreateGausSource_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_CreateTophatSource_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double BeamDiam,
    /* [in] */ long Nx,
    /* [in] */ long Ny,
    /* [in] */ double dx,
    /* [in] */ double dy,
    /* [in] */ double wavelength,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_CreateTophatSource_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ClearWavefront_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ClearWavefront_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMask_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double length_r,
    /* [in] */ double length_t,
    /* [in] */ double offset,
    /* [in] */ BSTR direction,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMask_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMaskMisalign_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double length_r,
    /* [in] */ double length_t,
    /* [in] */ double offset,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ BSTR direction,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMaskMisalign_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMaskGeneral_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double len_x,
    /* [in] */ double len_y,
    /* [in] */ double off_x,
    /* [in] */ double off_y,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMaskGeneral_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMaskPoly_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ SAFEARRAY __RPC_FAR * vertices,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMaskPoly_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMaskX_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double radius,
    /* [in] */ double length,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMaskX_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMaskXRounded_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double radius,
    /* [in] */ double length,
    /* [in] */ double corner_rad,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMaskXRounded_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMaskRotate_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double length_r,
    /* [in] */ double length_t,
    /* [in] */ double offset,
    /* [in] */ double angle,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ BSTR direction,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMaskRotate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ApplyMaskTilt_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double length_r,
    /* [in] */ double length_t,
    /* [in] */ double offset,
    /* [in] */ double rotangle,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ double tiltangle,
    /* [in] */ double tiltorientation,
    /* [in] */ BSTR direction,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ApplyMaskTilt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_Propagate_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double Distance,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_Propagate_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_PropagateExt_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double Distance,
    /* [in] */ double dxout,
    /* [in] */ double dyout,
    /* [in] */ long applyCurv,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_PropagateExt_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_PropagateDefDx_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double Distance,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_PropagateDefDx_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_CornerCube_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double ccsize,
    /* [in] */ BSTR shape,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ double spin,
    /* [in] */ SAFEARRAY __RPC_FAR * gapwidth,
    /* [in] */ SAFEARRAY __RPC_FAR * rotmatrix,
    /* [in] */ SAFEARRAY __RPC_FAR * dihedral,
    /* [in] */ SAFEARRAY __RPC_FAR * edgelength,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_CornerCube_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_focuslens_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double focuslens_f,
    /* [in] */ double focuslens_D,
    /* [in] */ double dxout,
    /* [in] */ double dyout,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_focuslens_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ThinLens_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double focallength,
    /* [in] */ double diameter,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ThinLens_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ReflectAsphere_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double rCurv,
    /* [in] */ double conicConst,
    /* [in] */ double diam,
    /* [in] */ double xDecenter,
    /* [in] */ double yDecenter,
    /* [in] */ double incidenceAngle,
    /* [in] */ double azimuth,
    /* [in] */ long applyCurv,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ReflectAsphere_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_PistonTiltFocus_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double piston,
    /* [in] */ double xtilt,
    /* [in] */ double ytilt,
    /* [in] */ double focus,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_PistonTiltFocus_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_PropagateToFocalPlane_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double lens_f,
    /* [in] */ double lens_D,
    /* [in] */ double dxout,
    /* [in] */ double dyout,
    /* [in] */ long applycurv,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_PropagateToFocalPlane_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ClipCirc_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double Diameter,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ double obscDiam,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ClipCirc_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ClipRect_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double len_x,
    /* [in] */ double len_y,
    /* [in] */ double off_x,
    /* [in] */ double off_y,
    /* [in] */ long obscFlag,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ClipRect_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ClipPolygon_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ SAFEARRAY __RPC_FAR * vertices,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ double angle,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ClipPolygon_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ClipBWindow_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ double winT,
    /* [in] */ double winA,
    /* [in] */ double xc,
    /* [in] */ double yc,
    /* [in] */ double rotAngle,
    /* [in] */ long BeamID,
    /* [retval][out] */ long __RPC_FAR *status);


void __RPC_STUB IDiffraction_ClipBWindow_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_AddThreadCommand_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long ThreadID,
    /* [in] */ long CmdID,
    /* [in] */ SAFEARRAY __RPC_FAR * rP,
    /* [in] */ SAFEARRAY __RPC_FAR * iP,
    /* [in] */ long pwAID,
    /* [in] */ long pwBID,
    /* [in] */ BSTR strP,
    /* [retval][out] */ long __RPC_FAR *flagStatus);


void __RPC_STUB IDiffraction_AddThreadCommand_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ExecuteThread_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long ThreadID,
    /* [retval][out] */ long __RPC_FAR *flagStatus);


void __RPC_STUB IDiffraction_ExecuteThread_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ThreadStatus_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long ThreadID,
    /* [retval][out] */ long __RPC_FAR *flagStatus);


void __RPC_STUB IDiffraction_ThreadStatus_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_WaitForThread_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long ThreadID,
    /* [in] */ long TimeOutMS,
    /* [retval][out] */ long __RPC_FAR *WaitStatus);


void __RPC_STUB IDiffraction_WaitForThread_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_WaitForMultipleThreads_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ SAFEARRAY __RPC_FAR * pThreadID,
    /* [in] */ long TimeOutMS,
    /* [in] */ long flagWaitAll,
    /* [retval][out] */ long __RPC_FAR *WaitStatus);


void __RPC_STUB IDiffraction_WaitForMultipleThreads_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [id] */ HRESULT STDMETHODCALLTYPE IDiffraction_ClearThread_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long ThreadID,
    /* [retval][out] */ long __RPC_FAR *flagStatus);


void __RPC_STUB IDiffraction_ClearThread_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [propget][id] */ HRESULT STDMETHODCALLTYPE IDiffraction_get_ThreadOPD_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long ThreadID,
    /* [in] */ long OPDID,
    /* [retval][out] */ double __RPC_FAR *VALUE);


void __RPC_STUB IDiffraction_get_ThreadOPD_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);


/* [propget][id] */ HRESULT STDMETHODCALLTYPE IDiffraction_get_ThreadPOW_Proxy( 
    IDiffraction __RPC_FAR * This,
    /* [in] */ long ThreadID,
    /* [in] */ long POWID,
    /* [retval][out] */ double __RPC_FAR *VALUE);


void __RPC_STUB IDiffraction_get_ThreadPOW_Stub(
    IRpcStubBuffer *This,
    IRpcChannelBuffer *_pRpcChannelBuffer,
    PRPC_MESSAGE _pRpcMessage,
    DWORD *_pdwStubPhase);



#endif 	/* __IDiffraction_INTERFACE_DEFINED__ */



#ifndef __DiffractionLib_LIBRARY_DEFINED__
#define __DiffractionLib_LIBRARY_DEFINED__

/* library DiffractionLib */
/* [helpstring][version][uuid] */ 


EXTERN_C const IID LIBID_DiffractionLib;

EXTERN_C const CLSID CLSID_Diffraction;

#ifdef __cplusplus

class DECLSPEC_UUID("E4BF1897-350C-488F-B246-E14BAF2A0756")
Diffraction;
#endif
#endif /* __DiffractionLib_LIBRARY_DEFINED__ */

/* Additional Prototypes for ALL interfaces */

unsigned long             __RPC_USER  BSTR_UserSize(     unsigned long __RPC_FAR *, unsigned long            , BSTR __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  BSTR_UserMarshal(  unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, BSTR __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  BSTR_UserUnmarshal(unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, BSTR __RPC_FAR * ); 
void                      __RPC_USER  BSTR_UserFree(     unsigned long __RPC_FAR *, BSTR __RPC_FAR * ); 

unsigned long             __RPC_USER  LPSAFEARRAY_UserSize(     unsigned long __RPC_FAR *, unsigned long            , LPSAFEARRAY __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  LPSAFEARRAY_UserMarshal(  unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, LPSAFEARRAY __RPC_FAR * ); 
unsigned char __RPC_FAR * __RPC_USER  LPSAFEARRAY_UserUnmarshal(unsigned long __RPC_FAR *, unsigned char __RPC_FAR *, LPSAFEARRAY __RPC_FAR * ); 
void                      __RPC_USER  LPSAFEARRAY_UserFree(     unsigned long __RPC_FAR *, LPSAFEARRAY __RPC_FAR * ); 

/* end of Additional Prototypes */

#ifdef __cplusplus
}
#endif

#endif
