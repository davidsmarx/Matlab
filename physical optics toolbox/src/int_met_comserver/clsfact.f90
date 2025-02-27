!
!  clsfact.f90 - This file contains methods of the IClassFactory 
!  interface that is used to create instances of the user's classes
!
!  Generated by the Visual Fortran COM Server Wizard on
!  08/08/05 at 17:24:41.
!
!   DO NOT EDIT THIS FILE!
!
!  This file is re-generated every time the object hierarchy is changed.
!

function IClassFactory_QueryInterface (pData, riid, ppv) result (r)
!DEC$ ATTRIBUTES STDCALL :: IClassFactory_QueryInterface
    use IClassFactory_Types
    use dfwinty
    use dfcom
    use kernel32
    implicit none

    type (IClassFactory_Data) pData
	!dec$ attributes reference :: pData
    type(GUID) riid
	!dec$ attributes reference :: riid
    integer(INT_PTR_KIND()) ppv
	!dec$ attributes reference :: ppv
    integer(LONG) r

    integer(LONG) i

    r = S_OK

    !  Ensure that they are requesting 
    !  the IClassFactory or IUnknown interface
    if ((.not. COMIsEqualGUID(riid, IID_IClassFactory)) .AND. &
        (.not. COMIsEqualGUID(riid, IID_IUnknown))) then
        r = E_UNEXPECTED 
        return
    end if

    i = InterlockedIncrement (loc(pData % refCount)) 
    ppv = loc(pData)

end function 

function IClassFactory_AddRef (pData) result (r)
!DEC$ ATTRIBUTES STDCALL :: IClassFactory_AddRef
    use IClassFactory_Types
    use dfwinty
    use kernel32
    implicit none

    type (IClassFactory_Data) pData
	!dec$ attributes reference :: pData
    integer r

    integer(LONG) i

    i = InterlockedIncrement (loc(pData % refCount)) 
    r = pData % refCount

end function 

function IClassFactory_Release (pData) result (r)
!DEC$ ATTRIBUTES STDCALL :: IClassFactory_Release
    use IClassFactory_Types
    use dfwinty
    use int_met_comserver_global
    use kernel32
    implicit none

    type (IClassFactory_Data), target :: pData
	!dec$ attributes reference :: pData
    integer r

    type (IClassFactory_Data),  pointer :: pCFData
    integer status

    call EnterCriticalSection(loc(gGlobalCriticalSection))
    pData % refCount = pData % refCount - 1
    r = pData % refCount
    if (pData % refCount == 0) then
        !  Time to delete ourself....
        deallocate (pData % pVtbl)
        pCFData => pData
        deallocate (pCFData)
        !status = ServerUnlock()
    end if
    call LeaveCriticalSection(loc(gGlobalCriticalSection))

end function 

function IClassFactory_LockServer (pData, bLock) result (r)
!DEC$ ATTRIBUTES STDCALL :: IClassFactory_LockServer
    use IClassFactory_Types
    use dfwinty
    use int_met_comserver_global
    implicit none

    type (IClassFactory_Data) pData
	!dec$ attributes reference :: pData
    logical bLock
	!dec$ attributes value :: bLock
    integer(LONG) r

    integer status

    r = S_OK

    if (bLock) then
        status = ServerLock()
    else
        status = ServerUnlock()
    end if

end function 

!  Per class
function IClassFactory_CreateDiffractionInstance (pData, pUnk, riid, ppv) result (r)
!DEC$ ATTRIBUTES STDCALL :: IClassFactory_CreateDiffractionInstance
    use IClassFactory_Types
    use dfwinty
    use dfcom
    use int_met_comserver_global
    use Diffraction_Types
    use IDiffraction_Methods
    implicit none

    type (IClassFactory_Data) pData
	!dec$ attributes reference :: pData
    integer(INT_PTR_KIND()) pUnk
	!dec$ attributes value :: pUnk
    type(GUID) riid
	!dec$ attributes reference :: riid
    integer(INT_PTR_KIND()) ppv
	!dec$ attributes reference :: ppv
    integer(LONG) r

    integer status
    integer(INT_PTR_KIND()) ptr
    type (Diffraction_Data),  pointer :: pDiffractionData

    r = S_OK
    ppv = NULL

    !  Ensure that they are requesting a supported interface
    if ((.not. COMIsEqualGUID(riid, IID_IUnknown)) &
        !  Per interface
         .AND. (.not. COMIsEqualGUID(riid, IID_IDiffraction)) &
         .AND. (.not. COMIsEqualGUID(riid, IID_IDispatch)) &
        ) then
        r = E_UNEXPECTED 
        return
    end if

	!  Class does not support aggregation
	if (pUnk /= NULL) then
        r = CLASS_E_NOAGGREGATION
        return
	endif

    !  Allocate instance data structure and Vtbl
    allocate (pDiffractionData)

    !  Per interface
    allocate (pDiffractionData % IDiffraction_InternalData % pVtbl)
    pDiffractionData % IDiffraction_InternalData % pVtbl % QueryInterface = loc(Diffraction_QueryInterface)
    pDiffractionData % IDiffraction_InternalData % pVtbl % AddRef = loc(Diffraction_AddRef)
    pDiffractionData % IDiffraction_InternalData % pVtbl % Release = loc(Diffraction_Release)
    pDiffractionData % IDiffraction_InternalData % pVtbl % GetTypeInfoCount = loc(Diffraction_GetTypeInfoCount)
    pDiffractionData % IDiffraction_InternalData % pVtbl % GetTypeInfo = loc(Diffraction_GetTypeInfo)
    pDiffractionData % IDiffraction_InternalData % pVtbl % GetIDsOfNames = loc(Diffraction_GetIDsOfNames)
    pDiffractionData % IDiffraction_InternalData % pVtbl % Invoke = loc(Diffraction_Invoke)
    pDiffractionData % IDiffraction_InternalData % pVtbl % x_vector_get = loc($IDiffraction_x_vector_get)
    pDiffractionData % IDiffraction_InternalData % pVtbl % y_vector_get = loc($IDiffraction_y_vector_get)
    pDiffractionData % IDiffraction_InternalData % pVtbl % wavefront_put = loc($IDiffraction_wavefront_put)
    pDiffractionData % IDiffraction_InternalData % pVtbl % wavefront_get = loc($IDiffraction_wavefront_get)
    pDiffractionData % IDiffraction_InternalData % pVtbl % CopyWavefront = loc($IDiffraction_CopyWavefront)
    pDiffractionData % IDiffraction_InternalData % pVtbl % CopyWavefrontAmp = loc($IDiffraction_CopyWavefrontAmp)
    pDiffractionData % IDiffraction_InternalData % pVtbl % SumWavefront = loc($IDiffraction_SumWavefront)
    pDiffractionData % IDiffraction_InternalData % pVtbl % GetWavefrontParms = loc($IDiffraction_GetWavefrontParms)
    pDiffractionData % IDiffraction_InternalData % pVtbl % WriteWavefrontUNF = loc($IDiffraction_WriteWavefrontUNF)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ReadWavefrontUNF = loc($IDiffraction_ReadWavefrontUNF)
    pDiffractionData % IDiffraction_InternalData % pVtbl % WavefrontPower = loc($IDiffraction_WavefrontPower)
    pDiffractionData % IDiffraction_InternalData % pVtbl % CreateGausSource = loc($IDiffraction_CreateGausSource)
    pDiffractionData % IDiffraction_InternalData % pVtbl % CreateTophatSource = loc($IDiffraction_CreateTophatSource)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ClearWavefront = loc($IDiffraction_ClearWavefront)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMask = loc($IDiffraction_ApplyMask)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMaskMisalign = loc($IDiffraction_ApplyMaskMisalign)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMaskGeneral = loc($IDiffraction_ApplyMaskGeneral)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMaskPoly = loc($IDiffraction_ApplyMaskPoly)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMaskX = loc($IDiffraction_ApplyMaskX)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMaskXRounded = loc($IDiffraction_ApplyMaskXRounded)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMaskRotate = loc($IDiffraction_ApplyMaskRotate)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ApplyMaskTilt = loc($IDiffraction_ApplyMaskTilt)
    pDiffractionData % IDiffraction_InternalData % pVtbl % Propagate = loc($IDiffraction_Propagate)
    pDiffractionData % IDiffraction_InternalData % pVtbl % PropagateExt = loc($IDiffraction_PropagateExt)
    pDiffractionData % IDiffraction_InternalData % pVtbl % PropagateDefDx = loc($IDiffraction_PropagateDefDx)
    pDiffractionData % IDiffraction_InternalData % pVtbl % CornerCube = loc($IDiffraction_CornerCube)
    pDiffractionData % IDiffraction_InternalData % pVtbl % focuslens = loc($IDiffraction_focuslens)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ThinLens = loc($IDiffraction_ThinLens)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ReflectAsphere = loc($IDiffraction_ReflectAsphere)
    pDiffractionData % IDiffraction_InternalData % pVtbl % PistonTiltFocus = loc($IDiffraction_PistonTiltFocus)
    pDiffractionData % IDiffraction_InternalData % pVtbl % PropagateToFocalPlane = loc($IDiffraction_PropagateToFocalPlane)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ClipCirc = loc($IDiffraction_ClipCirc)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ClipRect = loc($IDiffraction_ClipRect)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ClipPolygon = loc($IDiffraction_ClipPolygon)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ClipBWindow = loc($IDiffraction_ClipBWindow)
    pDiffractionData % IDiffraction_InternalData % pVtbl % AddThreadCommand = loc($IDiffraction_AddThreadCommand)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ExecuteThread = loc($IDiffraction_ExecuteThread)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ThreadStatus = loc($IDiffraction_ThreadStatus)
    pDiffractionData % IDiffraction_InternalData % pVtbl % WaitForThread = loc($IDiffraction_WaitForThread)
    pDiffractionData % IDiffraction_InternalData % pVtbl % WaitForMultipleThreads = loc($IDiffraction_WaitForMultipleThreads)
    pDiffractionData % IDiffraction_InternalData % pVtbl % ClearThread = loc($IDiffraction_ClearThread)
    pDiffractionData % IDiffraction_InternalData % pVtbl % get_ThreadOPD = loc($IDiffraction_get_ThreadOPD)
    pDiffractionData % IDiffraction_InternalData % pVtbl % get_ThreadPOW = loc($IDiffraction_get_ThreadPOW)
    pDiffractionData % IDiffraction_InternalData % pInternalData => pDiffractionData % InternalData
    pDiffractionData % InternalData % pStart => pDiffractionData
    pDiffractionData % InternalData % refCount = 1

    !  Allocate the user-defined instance data structure
    allocate (pDiffractionData % InternalData % pInstanceData)
    !  Call the "constructor"...
    status = Diffraction_CONSTRUCTOR(pDiffractionData % InternalData % pInstanceData)
    if (status /= S_OK) then
        deallocate (pDiffractionData % InternalData % pInstanceData)
        deallocate (pDiffractionData % IDiffraction_InternalData % pVtbl)
        deallocate (pDiffractionData)
        r = status
        return
    endif

    !  Increment the server's active object count
    status = ServerLock()

    !  Return the correct interface
    if (COMIsEqualGUID(riid, IID_IUnknown)) then
        ppv = loc(pDiffractionData)
    !  Per interface
    else if (COMIsEqualGUID(riid, IID_IDiffraction)) then
        ppv = loc(pDiffractionData % IDiffraction_InternalData)
    else if (COMIsEqualGUID(riid, IID_IDispatch)) then
        ppv = loc(pDiffractionData % IDiffraction_InternalData)
    end if

end function 
