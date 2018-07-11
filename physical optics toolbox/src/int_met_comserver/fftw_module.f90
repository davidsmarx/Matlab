! MODULE to make transparent calls to fftw, using the same interface as for Benson's
! fft calls
MODULE fftw_module

PUBLIC

    !! fftw constants copied from fftw3.f
      INTEGER, PARAMETER :: FFTW_R2HC = 0
      INTEGER, PARAMETER :: FFTW_HC2R = 1
      INTEGER, PARAMETER :: FFTW_DHT = 2
      INTEGER, PARAMETER :: FFTW_REDFT00=3
      INTEGER, PARAMETER :: FFTW_REDFT01=4
      INTEGER, PARAMETER :: FFTW_REDFT10=5
      INTEGER, PARAMETER :: FFTW_REDFT11=6
      INTEGER, PARAMETER :: FFTW_RODFT00=7
      INTEGER, PARAMETER :: FFTW_RODFT01=8
      INTEGER, PARAMETER :: FFTW_RODFT10=9
      INTEGER, PARAMETER :: FFTW_RODFT11=10
      INTEGER, PARAMETER :: FFTW_FORWARD = -1
      INTEGER, PARAMETER :: FFTW_BACKWARD=+1
      INTEGER, PARAMETER :: FFTW_MEASURE=0
      INTEGER, PARAMETER :: FFTW_DESTROY_INPUT=1
      INTEGER, PARAMETER :: FFTW_UNALIGNED=2
      INTEGER, PARAMETER :: FFTW_CONSERVE_MEMORY=4
      INTEGER, PARAMETER :: FFTW_EXHAUSTIVE=8
      INTEGER, PARAMETER :: FFTW_PRESERVE_INPUT=16
      INTEGER, PARAMETER :: FFTW_PATIENT=32
      INTEGER, PARAMETER :: FFTW_ESTIMATE=64
      INTEGER, PARAMETER :: FFTW_ESTIMATE_PATIENT=128
      INTEGER, PARAMETER :: FFTW_BELIEVE_PCOST=256
      INTEGER, PARAMETER :: FFTW_DFT_R2HC_ICKY=512
      INTEGER, PARAMETER :: FFTW_NONTHREADED_ICKY=1024
      INTEGER, PARAMETER :: FFTW_NO_BUFFERING=2048
      INTEGER, PARAMETER :: FFTW_NO_INDIRECT_OP=4096
      INTEGER, PARAMETER :: FFTW_ALLOW_LARGE_GENERIC=8192
      INTEGER, PARAMETER :: FFTW_NO_RANK_SPLITS=16384
      INTEGER, PARAMETER :: FFTW_NO_VRANK_SPLITS=32768
      INTEGER, PARAMETER :: FFTW_NO_VRECURSE=65536
      INTEGER, PARAMETER :: FFTW_NO_SIMD=131072

CONTAINS
    
SUBROUTINE FFTW1D( a, inverse, fftwplan )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calls for a 1-d FFTW of the input array.
   ! Uses the FFTW library. The FFTW library is best used by storing a "plan"
   ! and repeatedly calling the plan to execute the FFT. On the first call to
   ! FFTW1D, or when the SIZE(a) changes, the plan will be created. Optionally,
   ! the plan can be supplied. Note that the plan already knows if the fft is
   ! forward or backward. Also, it is assumed that all transforms are in-place,
   ! and that if a plan is supplied, it is for an in-place transform.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE KINDS

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:),INTENT(IN OUT) :: a
   LOGICAL,      OPTIONAL,    INTENT(IN)     :: inverse
   INTEGER*8,    OPTIONAL,    INTENT(IN)     :: fftwplan

   ! Local variables
   INTEGER*8, SAVE  :: planforward = 0
   INTEGER*8, SAVE  :: planbackward = 0
   INTEGER          :: fftdir
   INTEGER,   SAVE  :: n = 0
   COMPLEX(prc), POINTER, DIMENSION(:) :: atmp => NULL()


   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! if a plan is supplied, execute and go
   IF ( PRESENT( fftwplan ) ) THEN
        call dfftw_execute_dft(fftwplan, a, a)

   ELSE
        ! check if a new plan needs to be created
        IF (n .NE. SIZE(a) .OR. planforward .eq. 0) THEN
            n = SIZE(a)
            print*, "create new plan, n = ", n
            ! allocate a dummy array to make the plan because making the plan overwrites
            ! the array
            ALLOCATE( atmp(n) )
            call dfftw_plan_dft_1d(planforward,n,atmp,atmp,FFTW_FORWARD,FFTW_MEASURE)
            call dfftw_plan_dft_1d(planbackward,n,atmp,atmp,FFTW_BACKWARD,FFTW_MEASURE)
            DEALLOCATE( atmp )
        ENDIF

        !print*, "planforward = ", planforward
        !print*, "planbackward = ", planbackward
        !print*, "a(1) = ", a(1)
            
        ! check if inverse requested, and execute appropriate plan
        IF ( PRESENT( inverse ) ) THEN
            IF ( inverse ) THEN
                call dfftw_execute_dft(planbackward, a, a)
            ELSE 
                call dfftw_execute_dft(planforward, a, a)
            ENDIF
        ! else default is forward transform
        ELSE
            call dfftw_execute_dft(planforward, a, a)
            
        ENDIF
    ENDIF
    
END SUBROUTINE FFTW1D

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
SUBROUTINE FFTW2D( a, inverse, fftwplan )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calls for a 2-d FFTW of the input array.
   ! Uses the FFTW library. The FFTW library is best used by storing a "plan"
   ! and repeatedly calling the plan to execute the FFT. On the first call to
   ! FFTW1D, or when the SIZE(a) changes, the plan will be created. Optionally,
   ! the plan can be supplied. Note that the plan already knows if the fft is
   ! forward or backward. Also, it is assumed that all transforms are in-place,
   ! and that if a plan is supplied, it is for an in-place transform.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE KINDS

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:),INTENT(IN OUT) :: a
   LOGICAL,      OPTIONAL,    INTENT(IN)       :: inverse
   INTEGER*8,    OPTIONAL,    INTENT(IN)       :: fftwplan
 
   ! Local variables
   INTEGER*8, SAVE  :: planforward = 0
   INTEGER*8, SAVE  :: planbackward = 0
   INTEGER          :: fftdir
   INTEGER,   SAVE  :: nx = 0
   INTEGER,   SAVE  :: ny = 0
   COMPLEX(prc), POINTER, DIMENSION(:,:) :: atmp => NULL()

   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! if a plan is supplied, execute and go
   IF ( PRESENT( fftwplan ) ) THEN
        call dfftw_execute_dft(fftwplan, a, a)

   ELSE
        ! check if a new plan needs to be created
        IF (ny .NE. SIZE(a,1) .OR. nx .NE. SIZE(a,2) .OR. planforward .EQ. 0 ) THEN
            ny = SIZE(a,1)  ! rows
            nx = SIZE(a,2)  ! cols
            print*, "create a new plan, ny = ", ny, " nx = ", nx

            ! allocate a dummy array to make the plan because making the plan overwrites
            ! the array
            ALLOCATE( atmp(ny,nx) )

            call dfftw_plan_dft_2d(planforward,nx,ny,atmp,atmp,FFTW_FORWARD,FFTW_MEASURE)
            call dfftw_plan_dft_2d(planbackward,nx,ny,atmp,atmp,FFTW_BACKWARD,FFTW_MEASURE)

            DEALLOCATE( atmp )
        ENDIF

        ! check if inverse requested, and execute appropriate plan
        IF ( PRESENT( inverse ) ) THEN
            IF ( inverse ) THEN
                call dfftw_execute_dft(planbackward, a, a)
            ELSE 
                call dfftw_execute_dft(planforward, a, a)
            ENDIF
        ! else default is forward transform
        ELSE
            call dfftw_execute_dft(planforward, a, a)
        ENDIF
    ENDIF
    
END SUBROUTINE FFTW2D

END MODULE fftw_module
