!
!  UDiffractionTY.f90 - This module contains user-defined class 
!  definitions and methods
!

module Diffraction_USE

use dfmt, only : RTL_CRITICAL_SECTION, InitializeCriticalSection

USE Kinds
USE SI_Units
USE Wavefronts
USE ThreadWrappers

INTEGER, PARAMETER :: MAX_NR_THREADS = 6
INTEGER, PARAMETER :: MAX_NR_WFRONTS = 8
INTEGER, PARAMETER :: MAX_NR_OUTPUTS = 8

!! error codes
INTEGER, PARAMETER :: INTMETCOMNOERROR = 0
INTEGER, PARAMETER :: INVALIDTHREADID  = -10
INTEGER, PARAMETER :: INVALIDWFRONTID  = -11
INTEGER, PARAMETER :: INVALIDOUTPUTID  = -12
INTEGER, PARAMETER :: INVALIDWVSIZE    = -13
INTEGER, PARAMETER :: INVALIDGAPWIDTH  = -14
INTEGER, PARAMETER :: INVALIDDIHEDRAL  = -15
INTEGER, PARAMETER :: INVALIDEDGELENGTH= -16
INTEGER, PARAMETER :: INVALIDROTMATRIX = -17

    type Diffraction_InstanceData
        sequence
        !  TODO:  Add fields and remove "dummy"

        !! array of wavefronts for use in the command threads
        type(wavefront), DIMENSION(0:MAX_NR_WFRONTS) :: Wbeams

        !! array of threads:
        type(ExecuteCommandsThreadParameters), DIMENSION(MAX_NR_THREADS) :: Tcmds

        !! Critical Section Lock because fftw plan creation is not thread-safe
        type(RTL_CRITICAL_SECTION)                :: FFTWPlanLock
        !! nor is read and write wavefronts thread safe
        type(RTL_CRITICAL_SECTION)                :: ReadWriteUNFLock

        !! single threaded routines:
        !TYPE(PropagateThreadParameters)           :: TParms
        REAL(pr),     ALLOCATABLE, DIMENSION(:)   :: x_vector
        REAL(pr),     ALLOCATABLE, DIMENSION(:)   :: y_vector
        COMPLEX(prc), ALLOCATABLE, DIMENSION(:,:) :: beam

    
        !! dummy 4 byte variable provides alignment to allow sequence attribute
        INTEGER*4, DIMENSION(3)  :: dummy1

    end type

  contains

    !
    !  Constructor
    !
    function Diffraction_CONSTRUCTOR( ObjectData ) result (hresult)

        use dfwinty

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        integer(LONG) hresult

        integer ii

        hresult = S_OK
        !  TODO:  Add field initialization code

        ! PropagateThreadParameters TParms
        !NULLIFY(ObjectData%TParms%wvfront)
        !ObjectData%TParms%flagComplete = 1 ! initialize to complete state

        ! initialize the fftwlock critical section and each wavefront's pointer
        ! to the lock
        call InitializeCriticalSection( loc(ObjectData%FFTWPlanLock) )
        call InitializeCriticalSection( loc(ObjectData%ReadWriteUNFLock) )
        do ii = 1, MAX_NR_WFRONTS
            call AssignFFTWPlanLock( ObjectData%Wbeams(ii), loc(ObjectData%FFTWPlanLock) )
            call AssignReadWriteUNFLock( ObjectData%Wbeams(ii), loc(ObjectData%ReadWriteUNFLock) )
        end do ! for each wavefront
        
        ! ExecuteCommandsThreadParameters Tcmds
        do ii = 1, MAX_NR_THREADS
            NULLIFY(ObjectData%Tcmds(ii)%firstCommand)
            ObjectData%Tcmds(ii)%flagStatus = STATUS_NOERROR
            ALLOCATE(ObjectData%Tcmds(ii)%outputOPD(MAX_NR_OUTPUTS), &
                     ObjectData%Tcmds(ii)%outputPOW(MAX_NR_OUTPUTS) )
        end do

    end function

    !
    !  Destructor
    !
    subroutine Diffraction_DESTRUCTOR( ObjectData )
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        !  TODO:  Add field cleanup code
        integer ii

        do ii = 1, MAX_NR_THREADS
            DEALLOCATE(ObjectData%Tcmds(ii)%outputOPD, &
                     ObjectData%Tcmds(ii)%outputPOW )
        end do


    end subroutine

end module



