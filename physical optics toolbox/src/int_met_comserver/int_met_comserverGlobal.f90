!
!  int_met_comserverGlobal.f90 - This file contains the global data and functions
!  for the COM Server
!
!  Generated by the Visual Fortran COM Server Wizard on
!  08/08/05 at 17:24:41.
!
!   DO NOT EDIT THIS FILE!
!
!  This file is re-generated every time the object hierarchy is changed.
!

!
!  M O D U L E S
!

!  Module containg DLL global data

module int_met_comserver_Global
    use dfwinty
    use kernel32
    use IClassFactory_Types

    integer(HANDLE) ghinst
    type (T_RTL_CRITICAL_SECTION) gGlobalCriticalSection
    !  Per Class
    integer(INT_PTR_KIND()) gDiffractionIClassFactory
    type (IClassFactory_Data), pointer :: gpDiffractionCFData
    integer(INT_PTR_KIND()) gDiffractionITypeInfo

	! Type Library GUID		
	TYPE (GUID), PARAMETER :: GUID_TypeLib_int_met_comserver = &
		GUID(#7FB676A4, #8C3D, #49A4, &
        CHAR('80'X)//CHAR('4D'X)//CHAR('07'X)//CHAR('35'X)// &
        CHAR('30'X)//CHAR('B6'X)//CHAR('B4'X)//CHAR('00'X))

  contains

    integer function ServerLock ()
        use ole32
        integer(LONG) r
        r = CoAddRefServerProcess()
        ServerLock = r
    end function

    integer function ServerUnlock ()
        use ole32
        use user32
        integer(LONG) r
        r = CoReleaseServerProcess()
        if (r == 0) then
            call PostQuitMessage(0)
        end if
        ServerUnlock = r
    end function

end module int_met_comserver_Global

