!     Last change:  RSB  10 Oct 2002    6:12 pm
MODULE Wavefronts

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!
!  modified to use FFTW for all of the FFT calls
!  D. Marx
!  10/12/04
!============================

   USE Kinds
   USE Error_Exit
   USE Discrete_Transforms

   IMPLICIT NONE

   PRIVATE

   PRIVATE :: TestWavefrontAllocation

   PUBLIC :: ASSIGNMENT ( = )
!   PUBLIC :: OPERATOR   ( + )
!   PUBLIC :: OPERATOR   ( - )
!   PUBLIC :: OPERATOR   ( * )
!   PUBLIC :: CONJG

   PUBLIC :: WavefrontSum
   PUBLIC :: WavefrontDiff

   PUBLIC :: GenGaussBeam
   PUBLIC :: GenTopHatBeam
   PUBLIC :: GenFiberModeBeam

   PUBLIC :: WavefrontWavelength
   PUBLIC :: WavefrontPower
   PUBLIC :: NormalizeToPower
   PUBLIC :: SetWavefrontAmplitude
   PUBLIC :: ModifyCurvature
   PUBLIC :: ApplyCurvature
   PUBLIC :: ApplyTilt
   PUBLIC :: ApplyPiston
   PUBLIC :: PistonTiltFocus
   PUBLIC :: ApplyOPD
   PUBLIC :: SetZernikeAberration
   PUBLIC :: ApplyZernike
   PUBLIC :: ApplyAberration

   PUBLIC :: Apodize
   PUBLIC :: ClipCirc
   PUBLIC :: ClipRect
   PUBLIC :: ClipPoly
   PUBLIC :: SpiderObsc

   PUBLIC :: ClearWavefront
   PUBLIC :: GetWavefrontSampling
   PUBLIC :: GetWavefrontIntensity
   PUBLIC :: GetWavefrontAmp
   PUBLIC :: GetWavefrontOPD
   PUBLIC :: GetWavefrontSection
   PUBLIC :: AddWavefrontSection

   PUBLIC :: Propagate
   PUBLIC :: InvertWavefrontX
   PUBLIC :: InvertWavefrontY
   PUBLIC :: InvertWavefront
   PUBLIC :: ReflectFromCornerCube
   PUBLIC :: CoupleToFiber
   PUBLIC :: Heterodyne
   PUBLIC :: AverageOPD
   PUBLIC :: SaveWavefront
   PUBLIC :: RetrieveWavefront
   PUBLIC :: WriteWavefrontUNF

   ! DEFINITION OF WAVEFRONT OBJECT............................................
   ! Components initialized, in accord with Fortran 95
   TYPE, PUBLIC :: wavefront
      PRIVATE
      REAL(pr)                              :: dx = 0.0_pr     ! Sample spacing
      REAL(pr)                              :: dy = 0.0_pr
      REAL(pr)                              :: wavelength = 0.0_pr
      REAL(pr)                              :: curvature  = 0.0_pr
      COMPLEX(prc), POINTER, DIMENSION(:,:) :: amp => NULL()	! Field amplitude
      REAL(pr),     POINTER, DIMENSION(:)   :: x   => NULL()   ! Sample x-coords
      REAL(pr),     POINTER, DIMENSION(:)   :: y   => NULL()   ! Sample y-coords
      INTEGER*8                             :: fftw_plan2d_for, fftw_plan2d_rev
      TYPE(CztFFTWplans)                    :: cztplan
   END TYPE wavefront
   !...........................................................................

	! DEFINITION OF ZERNIKE ABERRATION OBJECT....................................
	INTEGER, PUBLIC, PARAMETER :: MAX_ZERNIKE_MODES = 36
	TYPE, PUBLIC :: zernikeAberration
		PRIVATE
		REAL(pr)                               :: refDiam    = 0.0_pr
		REAL(pr)                               :: xo         = 0.0_pr
		REAL(pr)                               :: yo         = 0.0_pr
		INTEGER,  DIMENSION(MAX_ZERNIKE_MODES) :: modeNum    = 0
		REAL(pr), DIMENSION(MAX_ZERNIKE_MODES) :: coeffValue = 0.0_pr
		LOGICAL                                :: isRMS      = .FALSE.
		INTEGER                                :: maxModes   = 0
	END TYPE zernikeAberration
	!............................................................................

   INTERFACE ASSIGNMENT ( = )
      MODULE PROCEDURE ReplaceWavefront
   END INTERFACE

!   INTERFACE OPERATOR ( + )
!      MODULE PROCEDURE AddWavefronts
!   END INTERFACE

!   INTERFACE OPERATOR ( - )
!      MODULE PROCEDURE SubtractWavefronts
!   END INTERFACE

!   INTERFACE OPERATOR ( * )
!      MODULE PROCEDURE MultiplyWavefronts
!   END INTERFACE

!   INTERFACE CONJG
!      MODULE PROCEDURE WavefrontConjugate
!   END INTERFACE

   INTERFACE GetWavefrontSection
      MODULE PROCEDURE GetWavefrontSection_XY, GetWavefrontSection_JK
   END INTERFACE

   INTERFACE ApplyZernike
      MODULE PROCEDURE ApplyOneZernike, ApplyZernikeSum
   END INTERFACE

	INTERFACE SetZernikeAberration
		MODULE PROCEDURE SetOneZernikeAberr, SetMultipleZernikeAberr
	END INTERFACE

   REAL(pr), PARAMETER :: GRID_CRIT       = 1.0D-6
   REAL(pr), PARAMETER :: WAVELENGTH_CRIT = 1.0D-9

CONTAINS


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ReplaceWavefront( beam, oldBeam )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Replaces beam with another wavefront (i.e., oldBeam): defines "=" for
   ! two wavefronts
   !     16 Oct  '00
   !     30 Nov  '01    Use DO loops for assigning values
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   TYPE(wavefront), INTENT(IN)     :: oldBeam

   INTRINSIC SIZE

   ! Local variables
   INTEGER :: Nxo, Nyo
   INTEGER :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! If possible, use the storage already allocated for beam; otherwise
   ! allocate new storage.

   Nxo = SIZE( oldBeam%x )
   Nyo = SIZE( oldBeam%y )

   CALL TestWavefrontAllocation( beam, Nxo, Nyo )

   beam%dx = oldBeam%dx
   beam%dy = oldBeam%dy
   beam%wavelength = oldBeam%wavelength
   beam%curvature  = oldBeam%curvature

   DO k = 1, Nxo
   	beam%x(k) = oldBeam%x(k)
   END DO
   DO j = 1,Nyo
   	beam%y(j) = oldBeam%y(j)
   END DO
   DO k = 1, Nxo
   	DO j = 1, Nyo
   		beam%amp(j,k) = oldBeam%amp(j,k)
   	END DO
   END DO

   END SUBROUTINE ReplaceWavefront



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE TestWavefrontAllocation( beam, Nx, Ny )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Allocates storage for wavefront parameters, if not already allocated
   ! Also creates fftw plans associated with the 2-d and 1-d complex arrays
   ! 
   !      6 Nov  '00
   !      2 Aug  '01    Simplified logic
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE fftw_module

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   INTEGER,         INTENT(IN)     :: Nx, Ny

   INTRINSIC ASSOCIATED, SIZE

   ! Local variables
   INTEGER            :: flag
   INTEGER, PARAMETER :: ok = 0
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ASSOCIATED( beam%amp ) ) THEN
      IF ( SIZE( beam%x ) == Nx  .AND.  SIZE( beam%y ) == Ny ) THEN
      	RETURN
      ELSE
        ! deallocate and clear plans
	   	DEALLOCATE( beam%x, beam%y, beam%amp, STAT=flag )
        call DFFTW_DESTROY_PLAN( beam%fftw_plan2d_for )
        call DFFTW_DESTROY_PLAN( beam%fftw_plan2d_rev )

        call DeAllocateCztFFTWplans( beam%cztplan )

        ! check success and clear other variables
   		IF ( flag /= ok )  &
      		CALL ErrorExit( "TestWavefrontAllocation: deallocation error" )
      	beam%dx = 0.0_pr
      	beam%dy = 0.0_pr
      	beam%wavelength = 0.0_pr
      	beam%curvature  = 0.0_pr
      END IF
   END IF

   ALLOCATE( beam%x(Nx), beam%y(Ny), beam%amp(Ny,Nx), STAT=flag )
   IF ( flag /= ok )  &
      CALL ErrorExit( "TestWavefrontAllocation: allocation error" )

   call dfftw_plan_dft_2d(beam%fftw_plan2d_for,nx,ny, &
        beam%amp,beam%amp,FFTW_FORWARD,FFTW_MEASURE)
   call dfftw_plan_dft_2d(beam%fftw_plan2d_rev,nx,ny, &
        beam%amp,beam%amp,FFTW_BACKWARD,FFTW_MEASURE)

   call AllocateCztFFTWplans( beam%cztplan, Nx, Ny )


   END SUBROUTINE TestWavefrontAllocation



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION AddWavefronts( beam1, beam2 ) RESULT( wavefrontSum )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Adds the complex values of two input wavefronts: defines "+"
   !
   !     CAUSES MEMORY LEAKS...USE SUBROUTINE WAVEFRONTSUM, INSTEAD
   !
   !     19 Oct  '00
   !     30 Mar  '01    Tests allocation of wavefrontSum
   !      2 Aug  '01    Tests wavelength difference; applies curvature
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines, ONLY : QuadraticPhase

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam1
   TYPE(wavefront), INTENT(IN) :: beam2
   TYPE(wavefront)             :: wavefrontSum

   INTRINSIC SIZE, ABS

   ! Local variables
   INTEGER  :: Ny, Nx
   REAL(pr) :: deltaLambda
   REAL(pr) :: piston, tilt, curvCoeff
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE( beam1%amp, 1 )
   Nx = SIZE( beam1%amp, 2 )

   deltaLambda = beam1%wavelength - beam2%wavelength

   IF ( Ny /= SIZE( beam2%amp, 1 ) .OR.  &
        Nx /= SIZE( beam2%amp, 2 ) ) THEN
      CALL ErrorExit( "AddWavefronts: unequal number of samples" )
   ELSE IF ( ABS( beam1%dx - beam2%dx ) > beam1%dx * GRID_CRIT  .OR. &
             ABS( beam1%dy - beam2%dy ) > beam1%dy * GRID_CRIT ) THEN
      CALL ErrorExit( "AddWavefronts: unequal sample spacing" )
   ELSE IF ( ABS( deltaLambda ) > beam1%wavelength*WAVELENGTH_CRIT ) THEN
      CALL ErrorExit( "AddWavefronts: unequal wavelengths" )
   ELSE
      CALL TestWavefrontAllocation( wavefrontSum, Nx, Ny )
      wavefrontSum%wavelength = beam2%wavelength
      wavefrontSum%dx  = beam2%dx
      wavefrontSum%dy  = beam2%dy
      wavefrontSum%x   = beam2%x
      wavefrontSum%y   = beam2%y
      wavefrontSum%amp = beam2%amp
      ! Apply curvature difference to beam2 before summing amplitudes
      piston = 0.0_pr
      tilt   = 0.0_pr
      curvCoeff = PI*(beam2%curvature - beam1%curvature)/beam2%wavelength
      CALL QuadraticPhase( wavefrontSum%amp, beam2%x, beam2%y, &
      							piston, tilt, tilt, curvCoeff )
      wavefrontSum%curvature = beam1%curvature
      wavefrontSum%amp       = beam1%amp + wavefrontSum%amp
   END IF

   END FUNCTION AddWavefronts


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE WavefrontSum( beam1, beam2, wfSum )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Adds the complex values of two input wavefronts: an alternative to the
   ! redefinition of the assignment operator (AddWavefronts), in order to
   ! avoid memory leaks
   !      2 Aug  '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines, ONLY : QuadraticPhase

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN)  :: beam1
   TYPE(wavefront), INTENT(IN)  :: beam2
   TYPE(wavefront), INTENT(IN OUT) :: wfSum

   INTRINSIC SIZE, ABS

   ! Local variables
   INTEGER  :: Ny, Nx
   REAL(pr) :: deltaLambda
   REAL(pr) :: piston, tilt, curvCoeff
   INTEGER  :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE( beam1%amp, 1 )
   Nx = SIZE( beam1%amp, 2 )

   deltaLambda = beam1%wavelength - beam2%wavelength

   IF ( Ny /= SIZE( beam2%amp, 1 ) .OR.  &
        Nx /= SIZE( beam2%amp, 2 ) ) THEN
      CALL ErrorExit( "AddWavefronts: unequal number of samples" )
   ELSE IF ( ABS( beam1%dx - beam2%dx ) > beam1%dx * GRID_CRIT  .OR. &
             ABS( beam1%dy - beam2%dy ) > beam1%dy * GRID_CRIT ) THEN
      CALL ErrorExit( "AddWavefronts: unequal sample spacing" )
   ELSE IF ( ABS( deltaLambda ) > beam1%wavelength*WAVELENGTH_CRIT ) THEN
      CALL ErrorExit( "AddWavefronts: unequal wavelengths" )
   ELSE


      !CALL ClearWavefront( wfSum )
      !ALLOCATE( wfSum%x(Nx), wfSum%y(Ny), wfSum%amp(Ny,Nx) )
      call TestWavefrontAllocation( wfSum, Nx, Ny )

      wfSum%wavelength = beam2%wavelength
      wfSum%dx  = beam2%dx
      wfSum%dy  = beam2%dy
      DO k = 1,Nx
      	wfSum%x(k) = beam2%x(k)
      END DO
      DO j = 1,Ny
      	wfSum%y(j) = beam2%y(j)
      END DO
      DO k = 1,Nx
      	DO j = 1,Ny
      		wfSum%amp(j,k) = beam2%amp(j,k)
      	END DO
      END DO
      ! Apply curvature difference to beam2 before summing amplitudes
      piston = 0.0_pr
      tilt   = 0.0_pr
      curvCoeff = PI*(beam2%curvature - beam1%curvature)/beam2%wavelength
      CALL QuadraticPhase( wfSum%amp, beam2%x, beam2%y, &
      							piston, tilt, tilt, curvCoeff )
      wfSum%curvature = beam1%curvature
      DO k = 1,Nx
      	DO j = 1,Ny
      		wfSum%amp(j,k) = beam1%amp(j,k) + wfSum%amp(j,k)
      	END DO
      END DO
   END IF

   END SUBROUTINE WavefrontSum



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE WavefrontDiff( beam1, beam2, wfDiff )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Subtracts the complex values of two input wavefronts: an alternative to
   ! the redefinition of the assignment operator (SubtractWavefronts), in order
   ! to avoid memory leaks
   !      2 Aug  '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines, ONLY : QuadraticPhase

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN)  :: beam1
   TYPE(wavefront), INTENT(IN)  :: beam2
   TYPE(wavefront), INTENT(IN OUT) :: wfDiff

   INTRINSIC SIZE, ABS

   ! Local variables
   INTEGER  :: Ny, Nx
   REAL(pr) :: deltaLambda
   REAL(pr) :: piston, tilt, curvCoeff
   INTEGER  :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE( beam1%amp, 1 )
   Nx = SIZE( beam1%amp, 2 )

   deltaLambda = beam1%wavelength - beam2%wavelength

   IF ( Ny /= SIZE( beam2%amp, 1 ) .OR.  &
        Nx /= SIZE( beam2%amp, 2 ) ) THEN
      CALL ErrorExit( "WavefrontDiff: unequal number of samples" )
   ELSE IF ( ABS( beam1%dx - beam2%dx ) > beam1%dx * GRID_CRIT  .OR. &
             ABS( beam1%dy - beam2%dy ) > beam1%dy * GRID_CRIT ) THEN
      CALL ErrorExit( "WavefrontDiff: unequal sample spacing" )
   ELSE IF ( ABS( deltaLambda ) > beam1%wavelength*WAVELENGTH_CRIT ) THEN
      CALL ErrorExit( "WavefrontDiff: unequal wavelengths" )
   ELSE

      !CALL ClearWavefront( wfDiff )
      !ALLOCATE( wfDiff%x(Nx), wfDiff%y(Ny), wfDiff%amp(Ny,Nx) )
      call TestWavefrontAllocation( wfDiff, Nx, Ny )


      wfDiff%wavelength = beam2%wavelength
      wfDiff%dx  = beam2%dx
      wfDiff%dy  = beam2%dy
      DO k = 1,Nx
      	wfDiff%x(k) = beam2%x(k)
      END DO
      DO j = 1,Ny
      	wfDiff%y(j) = beam2%y(j)
      END DO
      DO k = 1,Nx
      	DO j = 1,Ny
      		wfDiff%amp(j,k) = beam2%amp(j,k)
      	END DO
      END DO
      ! Apply curvature difference to beam2 before subtracting amplitudes
      piston = 0.0_pr
      tilt   = 0.0_pr
      curvCoeff = PI*(beam2%curvature - beam1%curvature)/beam2%wavelength
      CALL QuadraticPhase( wfDiff%amp, beam2%x, beam2%y, &
      							piston, tilt, tilt, curvCoeff )
      wfDiff%curvature = beam1%curvature
      DO k = 1,Nx
      	DO j = 1,Ny
      		wfDiff%amp(j,k) = beam1%amp(j,k) - wfDiff%amp(j,k)
      	END DO
      END DO
   END IF

   END SUBROUTINE WavefrontDiff





  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION SubtractWavefronts( beam1, beam2 ) RESULT( wavefrontDiff )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Subtracts the complex values of two input wavefronts: defines "-"
   !
   !     PROBABLY CAUSES MEMORY LEAKS...USE SUBROUTINE WAVEFRONTDIFFERENCE
   !
   !     19 Oct  '00
   !     30 Mar  '01    Tests allocation of wavefrontDiff
   !      2 Aug  '01    Tests wavelength difference; applies curvature
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines, ONLY : QuadraticPhase

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam1
   TYPE(wavefront), INTENT(IN) :: beam2
   TYPE(wavefront)             :: wavefrontDiff

   INTRINSIC SIZE, ABS

   ! Local variables
   INTEGER  :: Ny, Nx
   REAL(pr) :: deltaLambda
   REAL(pr) :: piston, tilt, curvCoeff
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE( beam1%amp, 1 )
   Nx = SIZE( beam1%amp, 2 )

   deltaLambda = beam1%wavelength - beam2%wavelength

   IF ( Ny /= SIZE( beam2%amp, 1 )  .OR.  &
        Nx /= SIZE( beam2%amp, 2 ) ) THEN
      CALL ErrorExit( "SubtractWavefronts: unequal number of samples" )
   ELSE IF ( ABS( beam1%dx - beam2%dx ) > beam1%dx * GRID_CRIT  .OR. &
             ABS( beam1%dy - beam2%dy ) > beam1%dy * GRID_CRIT ) THEN
      CALL ErrorExit( "SubtractWavefronts: unequal sample spacing" )
   ELSE IF ( ABS( deltaLambda ) > beam1%wavelength*WAVELENGTH_CRIT ) THEN
      CALL ErrorExit( "SubtractWavefronts: unequal wavelengths" )
   ELSE
      wavefrontDiff = beam2
      ! Apply curvature difference to beam2 before subtracting amplitudes
      piston = 0.0_pr
      tilt   = 0.0_pr
      curvCoeff = PI*(beam2%curvature - beam1%curvature)/beam2%wavelength
      CALL QuadraticPhase( wavefrontDiff%amp, beam2%x, beam2%y, &
      							piston, tilt, tilt, curvCoeff )
      wavefrontDiff%curvature = beam1%curvature
      wavefrontDiff%amp       = beam1%amp - wavefrontDiff%amp
   END IF

   END FUNCTION SubtractWavefronts



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION MultiplyWavefronts( beam1, beam2 ) RESULT( wavefrontProd )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Multiplies the complex values of two input wavefronts, point-by-point:
   ! defines "*" for wavefronts
   !
   !     PROBABLY CAUSES MEMORY LEAKS...
   !
   !     19 Oct  '00
   !     19 Mar  '01    Corrected error message
   !     30 Mar  '01    Tests allocation of wavefrontProd
   !      2 Aug  '01    Tests wavelength difference; adds curvatures
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam1
   TYPE(wavefront), INTENT(IN) :: beam2
   TYPE(wavefront)             :: wavefrontProd

   INTRINSIC SIZE, ABS

   ! Local variables
   INTEGER  :: Ny, Nx
   REAL(pr) :: deltaLambda
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE( beam1%amp, 1 )
   Nx = SIZE( beam1%amp, 2 )

   deltaLambda = beam1%wavelength - beam2%wavelength

   IF ( Ny /= SIZE( beam2%amp, 1 )  .OR.  &
        Nx /= SIZE( beam2%amp, 2 ) ) THEN
      CALL ErrorExit( "MultiplyWavefronts: unequal number of samples" )
   ELSE IF ( ABS( beam1%dx - beam2%dx ) > beam1%dx * GRID_CRIT  .OR. &
             ABS( beam1%dy - beam2%dy ) > beam1%dy * GRID_CRIT ) THEN
      CALL ErrorExit( "MultiplyWavefronts: unequal sample spacing" )
   ELSE IF ( ABS( deltaLambda ) > beam1%wavelength*WAVELENGTH_CRIT ) THEN
      CALL ErrorExit( "MultiplyWavefronts: unequal wavelengths" )
   ELSE
      wavefrontProd = beam1
      wavefrontProd%amp       = wavefrontProd%amp * beam2%amp
      wavefrontProd%curvature = beam1%curvature + beam2%curvature
   END IF

   END FUNCTION MultiplyWavefronts



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION WavefrontConjugate( beam ) RESULT( conjugate )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns the complex conjugate of an input wavefront: defines CONJG() for
   ! wavefront argument
   !
   !     PROBABLY CAUSES MEMORY LEAKS...
   !
   !     19 Oct  '00
   !     30 Mar  '01    Tests allocation of conjugate (NOT NECESSARY...Aug '01)
   !      2 Aug  '01    Change sign of curvature
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   TYPE(wavefront)             :: conjugate

   INTRINSIC CONJG

   ! Local variables
   INTEGER :: Nx, Ny
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE( beam%amp, 1 )
   Nx = SIZE( beam%amp, 2 )

   conjugate = beam

   conjugate%amp       = CONJG( conjugate%amp )
   conjugate%curvature = -beam%curvature

   END FUNCTION WavefrontConjugate



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ClearWavefront( beam )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Free the storage used by a wavefront
   !     16 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam

   INTRINSIC ASSOCIATED
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ASSOCIATED( beam%amp ) ) THEN
      beam%dx = 0.0_pr
      beam%dy = 0.0_pr
      beam%wavelength = 0.0_pr
      beam%curvature  = 0.0_pr


      DEALLOCATE( beam%x, beam%y, beam%amp )
      
      call DFFTW_DESTROY_PLAN( beam%fftw_plan2d_for )
      call DFFTW_DESTROY_PLAN( beam%fftw_plan2d_rev )
      call DeAllocateCztFFTWplans( beam%cztplan )

   END IF

   END SUBROUTINE ClearWavefront


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GetWavefrontSampling( beam, Nx, Ny, dx, dy, xo, yo )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns the number of transverse samples in the specified wavefront,
   ! and their separation
   !     18 Oct  '00
   !     22 Apr  '02    All output quantities are optional; added xo and yo

   !  beam     A wavefront
   !  Nx, Ny   Number of samples in x- and y-direction
   !  dx, dy   Sample spacing
   !  xo, yo   Coordinate value of grid center point
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront),    INTENT(IN)  :: beam
   INTEGER,  OPTIONAL, INTENT(OUT) :: Nx, Ny
   REAL(pr), OPTIONAL, INTENT(OUT) :: dx, dy
   REAL(pr), OPTIONAL, INTENT(OUT) :: xo, yo

   INTRINSIC ASSOCIATED, SIZE, PRESENT
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ASSOCIATED( beam%amp ) ) THEN
      IF ( PRESENT( Nx ) ) Nx = SIZE( beam%amp, 2 )
      IF ( PRESENT( Ny ) ) Ny = SIZE( beam%amp, 1 )
      IF ( PRESENT( dx ) ) dx = beam%dx
      IF ( PRESENT( dy ) ) dy = beam%dy
      IF ( PRESENT( xo ) ) xo = beam%x(1+(Nx/2))
      IF ( PRESENT( yo ) ) yo = beam%y(1+(Ny/2))
   ELSE
      IF ( PRESENT( Nx ) ) Nx = 0
      IF ( PRESENT( Ny ) ) Ny = 0
      IF ( PRESENT( dx ) ) dx = 0.0_pr
      IF ( PRESENT( dy ) ) dy = 0.0_pr
      IF ( PRESENT( xo ) ) xo = 0.0_pr
      IF ( PRESENT( yo ) ) yo = 0.0_pr
   END IF


   END SUBROUTINE GetWavefrontSampling



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GenGaussBeam( Ugauss, Nx, Ny, dx, dy, lambda, Dbeam,  &
                            xc, yc, Dclip, power )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Create a Gaussian beam.
   ! Initial curvature is zero (i.e., flat wavefront)
   !  30 Oct  '00
   !  13 Dec  '00 Subroutine, rather than function
   !  10 Aug  '01 Made amplitude calculation a bit more tidy
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : CoordinateVector

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT)  :: Ugauss

   INTEGER,  INTENT(IN)           :: Nx, Ny           ! Number of samples
   REAL(pr), INTENT(IN)           :: dx, dy           ! Sample spacing
   REAL(pr), INTENT(IN)           :: lambda           ! Beam wavelength
   REAL(pr), INTENT(IN)           :: Dbeam            ! Twice beam spot size
   REAL(pr), INTENT(IN), OPTIONAL :: xc, yc           ! Beam center
   REAL(pr), INTENT(IN), OPTIONAL :: Dclip            ! Diameter for clipping
   REAL(pr), INTENT(IN), OPTIONAL :: power            ! Beam power

   INTRINSIC PRESENT, CMPLX, EXP

   !  Local variables
   INTEGER  :: j, k
   REAL(pr) :: spot, xbc, ybc
   REAL(pr) :: xsq,  ysq, realPart
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xbc = 0.0_pr
   ybc = 0.0_pr
   IF ( PRESENT( xc ) ) xbc = xc
   IF ( PRESENT( yc ) ) ybc = yc
   spot = ( 0.5_pr * Dbeam )**2

   !CALL ClearWavefront( Ugauss )
   !ALLOCATE( Ugauss%x(Nx), Ugauss%y(Ny), Ugauss%amp(Ny,Nx) )
   call TestWavefrontAllocation( Ugauss, Nx, Ny )

   Ugauss%dx = dx
   Ugauss%dy = dy
   Ugauss%x = CoordinateVector( Nx, dx )
   Ugauss%y = CoordinateVector( Ny, dy )

   DO k = 1,Nx
      xsq = (Ugauss%x(k)-xbc)**2
   	DO j = 1,Ny
      	ysq = (Ugauss%y(j)-ybc)**2
      	realPart = EXP( -( ( xsq + ysq )/spot ) )
      	Ugauss%amp(j,k) = CMPLX( realPart, 0.0_pr, KIND=prc )
   	END DO
   END DO

   Ugauss%wavelength = lambda
   Ugauss%curvature  = 0.0_pr

   IF ( PRESENT( Dclip ) ) CALL ClipCirc( Ugauss, Dclip )
   IF ( PRESENT( power ) ) CALL NormalizeToPower( Ugauss, power )

   END SUBROUTINE GenGaussBeam



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GenTopHatBeam( TopHat, Nx, Ny, dx, dy, lambda, Dbeam,  &
                             xc, yc, power )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Create a top-hat beam
   ! Initial curvature is zero (i.e., flat wavefront)
   !     30 Oct  '00
   !  	13 Dec  '00 Subroutine, rather than function
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : CoordinateVector

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: TopHat

   INTEGER,  INTENT(IN)           :: Nx, Ny           ! Number of samples
   REAL(pr), INTENT(IN)           :: dx, dy           ! Sample spacing
   REAL(pr), INTENT(IN)           :: lambda           ! Beam wavelength
   REAL(pr), INTENT(IN)           :: Dbeam            ! Beam diameter
   REAL(pr), INTENT(IN), OPTIONAL :: xc, yc           ! Beam center
   REAL(pr), INTENT(IN), OPTIONAL :: power            ! Beam power

   INTRINSIC PRESENT

   !  Local variables
   INTEGER  :: j, k
   COMPLEX(prc), PARAMETER :: ONE = ( 1.0_pr, 0.0_pr )
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   !CALL ClearWavefront( TopHat )
   !ALLOCATE( TopHat%x(Nx), TopHat%y(Ny), TopHat%amp(Ny,Nx) )
   call TestWavefrontAllocation( TopHat, Nx, Ny )

   TopHat%dx = dx
   TopHat%dy = dy
   TopHat%x = CoordinateVector( Nx, dx )
   TopHat%y = CoordinateVector( Ny, dy )

   TopHat%amp = ONE

   TopHat%wavelength = lambda
   TopHat%curvature  = 0.0_pr

   CALL ClipCirc( TopHat, Dbeam, xc, yc )

   IF ( PRESENT( power ) ) CALL NormalizeToPower( TopHat, power )

   END SUBROUTINE GenTopHatBeam



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GenFiberModeBeam( Ufiber, Nx, Ny, dx, dy, lambda,  &
										  coreDiam, numAper, xc, yc, power )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Create a beam based on the E-field distribution for the HE(11) mode of
   ! a step-index optical fiber.  Curvature is zero (i.e., flat wavefront).
   ! Default normalization (if "power" is not present) gives unit self-overlap.

   !   6 Nov  '00
   !  13 Dec  '00 Subroutine, rather than function
   !   2 Nov '01  Calls FiberMode subroutine (not function)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : CoordinateVector

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: Ufiber

   INTEGER,  INTENT(IN)           :: Nx, Ny           ! Number of samples
   REAL(pr), INTENT(IN)           :: dx, dy           ! Sample spacing
   REAL(pr), INTENT(IN)           :: lambda           ! Beam wavelength
   REAL(pr), INTENT(IN)           :: coreDiam         ! Diameter of fiber core
   REAL(pr), INTENT(IN)           :: numAper          ! Numerical aperture
   REAL(pr), INTENT(IN), OPTIONAL :: xc, yc           ! Location of core center
   REAL(pr), INTENT(IN), OPTIONAL :: power            ! Beam power

   INTRINSIC PRESENT

   !  Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   !CALL ClearWavefront( Ufiber )
   !ALLOCATE( Ufiber%x(Nx), Ufiber%y(Ny), Ufiber%amp(Ny,Nx) )
   call TestWavefrontAllocation( Ufiber, Nx, Ny )

   Ufiber%dx = dx
   Ufiber%dy = dy
   Ufiber%x = CoordinateVector( Nx, dx )
   Ufiber%y = CoordinateVector( Ny, dy )

   CALL FiberMode( Ufiber%x, Ufiber%y, lambda,           &
                   coreDiam, numAper, xc, yc, Ufiber%amp )
   Ufiber%wavelength = lambda
   Ufiber%curvature  = 0.0_pr

   IF ( PRESENT( power ) ) CALL NormalizeToPower( Ufiber, power )

   END SUBROUTINE GenFiberModeBeam



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FiberMode( x, y, lambda, coreDiam, numAper, xc, yc, modeAmp )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns the E-field amplitude for the HE(11) mode of a single-mode fiber.
   ! See: Peter K. Cheo, "Fiber Optics & Optoelectronics," 2nd ed.,
   !                     Prentice Hall, N.J. (1990), chapters 4 and 6.
   ! Curvature is zero (i.e., flat wavefront).
   ! Normalization is for unit self-overlap.

   !      6 Nov  '00
   !      2 Nov  '01    Changed to subroutine, to minimize stack storage
   !     30 Nov  '01    Use DO loops to normalize mode (save temp storage)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Math_Routines

   IMPLICIT NONE

   REAL(pr), DIMENSION(:), INTENT(IN) :: x         ! X-coordinate values
   REAL(pr), DIMENSION(:), INTENT(IN) :: y         ! Y-coordinate values
   REAL(pr), INTENT(IN)               :: lambda    ! Beam wavelength
   REAL(pr), INTENT(IN)               :: coreDiam  ! Diameter of fiber core
   REAL(pr), INTENT(IN)               :: numAper   ! Numerical aperture
   REAL(pr), INTENT(IN), OPTIONAL     :: xc, yc    ! Location of core center

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: modeAmp  ! HE11 mode

   INTRINSIC PRESENT, SIZE, SQRT, CMPLX

   !  Local variables
   REAL(pr) :: xbc, ybc
   REAL(pr) :: V, u, g, ur, gr
   REAL(pr) :: coreRadius
   REAL(pr) :: amp, beamSum, const, xsq, r
   INTEGER  :: j,   k
   INTEGER  :: Nx,  Ny
   REAL(pr) :: dx,  dy
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( modeAmp, 2 )
   Ny = SIZE( modeAmp, 1 )
   IF ( Nx /= SIZE( x )  .OR.  Ny /= SIZE( y ) )  &
      CALL ErrorExit( "FiberMode: array sizes do not match" )

   xbc = 0.0_pr
   ybc = 0.0_pr
   IF ( PRESENT( xc ) ) xbc = xc
   IF ( PRESENT( yc ) ) ybc = yc

   V = PI*coreDiam*numAper/lambda
   IF ( V > 2.405_pr ) CALL ErrorExit( "FiberMode: fiber is not single-mode" )

   u = ( 1.0_pr + sqrt2 )*V / ( 1.0_pr + SQRT( SQRT( 4.0_pr + V**4 ) ) )
   g = SQRT( V**2 - u**2 )

   coreRadius = 0.5_pr * coreDiam
   ur = u/coreRadius
   gr = g/coreRadius

   const = BesselJ0( u ) / BesselK0( g )

   beamSum = 0.0_pr

   DO k = 1,Nx
      xsq = ( x(k) - xbc )**2
   	DO j = 1,Ny
         r = SQRT( xsq + ( y(j) - ybc )**2 )
         IF ( r <= coreRadius ) THEN
            amp = BesselJ0( ur*r )
         ELSE
            amp = BesselK0( gr*r )
         END IF
         beamSum = beamSum + amp**2
         modeAmp(j,k) = CMPLX( amp, 0.0_pr, KIND=prc )
      END DO
   END DO

   dx = ( x(Nx) - x(1) )/REAL( Nx-1, KIND=pr )
   dy = ( y(Ny) - y(1) )/REAL( Ny-1, KIND=pr )

   DO k = 1, Nx
   	DO j = 1, Ny
   		modeAmp(j,k) = modeAmp(j,k) / (dx*dy*beamSum)
   	END DO
   END DO

   END SUBROUTINE FiberMode



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION CoupleToFiber( beam, coreDiam, numAper, xc, yc ) RESULT( coeff )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates overlap integral of beam and fiber mode; returns complex-valued
   ! coupling coefficient

   !      6 Nov  '00
   !      2 Nov  '01    Calls FiberMode subroutine (not function)
   !     30 Nov  '01    Use DO loops to save temp storage
   !     17 Dec  '01    Multiply (not divide) coeff by (dx*dy)--bug fix
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr),        INTENT(IN)     :: coreDiam
   REAL(pr),        INTENT(IN)     :: numAper
   REAL(pr),        INTENT(IN), OPTIONAL :: xc, yc

   COMPLEX(prc) :: coeff

   INTRINSIC SIZE, CONJG

   ! Local variables
   INTEGER :: j, k
   COMPLEX(prc), DIMENSION(SIZE(beam%y),SIZE(beam%x)) :: modeAmp
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   CALL FiberMode( beam%x, beam%y, beam%wavelength,   &
   					 coreDiam, numAper, xc, yc, modeAmp )

   coeff = ( 0.0_pr, 0.0_pr )
   DO k = 1, SIZE( beam%x )
   	DO j = 1, SIZE( beam%y )
   		coeff = coeff + beam%amp(j,k) * CONJG( modeAmp(j,k) )
   	END DO
   END DO

   coeff = coeff*(beam%dx*beam%dy)

   END FUNCTION CoupleToFiber



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION WavefrontWavelength( beam ) RESULT( lambda )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns value of wavelength for the input wavefront

   !      7 Aug '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   REAL(pr)                    :: lambda

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   lambda = beam%wavelength

   END FUNCTION WavefrontWavelength



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION WavefrontPower( beam ) RESULT( power )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates power of the input beam.
   !     27 May  '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   REAL(pr)                    :: power

   INTRINSIC SIZE, REAL, AIMAG

   ! Local variables
   INTEGER  :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   power = 0.0_pr
   DO k = 1,SIZE( beam%amp, 2 )
   	DO j = 1,SIZE( beam%amp, 1 )
      	power = power + REAL( beam%amp(j,k) )**2 + AIMAG( beam%amp(j,k) )**2
   	END DO
   END DO

   power = power * (beam%dx*beam%dy)

   END FUNCTION WavefrontPower



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE NormalizeToPower( beam, power )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Normalizes the input beam to the specified power
   !     16 Oct  '00
   !     27 May  '01    Calls WavefrontPower
   !     30 Nov  '01    Use DO loops to save temp storage
   !     17 Dec  '01    Corrected bug (j and k indexes for beam%amp in loop)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr),        INTENT(IN)     :: power

   INTRINSIC SQRT, EPSILON, SIZE

   ! Local variables
   REAL(pr) :: origPower, factor
   INTEGER  :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   origPower = WavefrontPower( beam )

   IF ( origPower > EPSILON( 0.0_pr ) ) THEN
      factor    = SQRT( power/origPower )
      DO k = 1, SIZE( beam%amp, 2 )
      	DO j = 1, SIZE( beam%amp, 1 )
      		beam%amp(j,k)  = beam%amp(j,k) * factor
      	END DO
      END DO
   ELSE
      CALL ErrorExit( "NormalizeToPower: zero-power beam?" )
   END IF

   END SUBROUTINE NormalizeToPower



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SetWavefrontAmplitude( beam, amplitudeValue )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Sets all values of amplitude to the specified (complex) input value.

   !      8 Apr '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   COMPLEX(prc),    INTENT(IN)     :: amplitudeValue

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   beam%amp = amplitudeValue

   END SUBROUTINE SetWavefrontAmplitude



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ModifyCurvature( beam, deltaCurv )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Changes the value of wavefront curvature for the input wavefront.
   ! For example, a simple lens of focal length f causes a curvature
   ! change of -1/f (although "ThinLens" is a better way to represent a lens).
   !     16 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr),        INTENT(IN)     :: deltaCurv
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   beam%curvature = beam%curvature + deltaCurv

   END SUBROUTINE ModifyCurvature



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ApplyCurvature( beam, curvatureValue, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Explicitly applies the beam's curvature to the array of field values.
   ! However, if curvatureValue is present, the input curvature value
   ! is applied to the wavefront; in this case, optionally, the center of
   ! curvature may be specified by (xc,yc).

   !     20 Nov  '00
   !     16 July '01    Revised to call QuadraticPhase
   !      4 Jan  '02    Revised to not require x() and y() work arrays
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines

   IMPLICIT NONE

   TYPE(wavefront),    INTENT(IN OUT) :: beam
   REAL(pr), OPTIONAL, INTENT(IN)     :: curvatureValue
   REAL(pr), OPTIONAL, INTENT(IN)     :: xc, yc

   INTRINSIC PRESENT

   ! Local variables
   REAL(pr) :: curvCoeff, xTilt, yTilt, piston
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xTilt  = 0.0_pr
   yTilt  = 0.0_pr
   piston = 0.0_pr

   IF ( PRESENT( curvatureValue ) )THEN
      ! Apply input curvature value to complex field
      curvCoeff = PI*curvatureValue/beam%wavelength
      IF ( PRESENT( xc ) ) THEN
   	   xTilt  = -2.0_pr*xc*curvCoeff
         piston = curvCoeff * xc**2
      END IF
      IF ( PRESENT( yc ) ) THEN
   	   yTilt  = -2.0_pr*yc*curvCoeff
         piston = curvCoeff * yc**2 + piston
      END IF
   ELSE
      ! Apply beam's own curvature to complex field
      curvCoeff = PI*beam%curvature/beam%wavelength
      beam%curvature = 0.0_pr
   END IF

   CALL QuadraticPhase( beam%amp, beam%x, beam%y,       &
   							piston, xTilt, yTilt, curvCoeff )

   END SUBROUTINE ApplyCurvature



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ApplyTilt( beam, xTilt, yTilt, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Explicitly applies tilt to the array of field values
   ! -- xTilt and yTilt are angles in radians --
   !     19 Oct  '00
   !     16 July '01    Revised to call QuadraticPhase
   !     12 Dec  '01    Revised to not require x() and y() work arrays
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines

   IMPLICIT NONE

   TYPE(wavefront),    INTENT(IN OUT)       :: beam
   REAL(pr),           INTENT(IN)           :: xTilt, yTilt
   REAL(pr),           INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC PRESENT

   ! Local variables
   REAL(pr) :: xTiltCoeff, yTiltCoeff, piston, curv
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xTiltCoeff = TWOPI*xTilt/beam%wavelength
   yTiltCoeff = TWOPI*yTilt/beam%wavelength

   piston = 0.0_pr
   curv   = 0.0_pr

   IF ( PRESENT( xc ) ) piston = -xc*xTiltCoeff
   IF ( PRESENT( yc ) ) piston = -yc*yTiltCoeff + piston

   CALL QuadraticPhase( beam%amp, beam%x, beam%y,            &
   							piston, xTiltCoeff, yTiltCoeff, curv )

   END SUBROUTINE ApplyTilt



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ApplyPiston( beam, piston )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Explicitly applies piston to the array of field values
   ! Piston input is OPD in meters.
   !     10 Nov  '00
   !     30 Nov  '01    Use DO loops to save temp storage
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr),        INTENT(IN)     :: piston

   INTRINSIC SIZE, EXP, CMPLX

   ! Local variables
   REAL(pr)     :: phase
   COMPLEX(prc) :: const
   INTEGER      :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   phase = TWOPI*piston/beam%wavelength
   const = EXP( CMPLX( 0.0_pr, phase, KIND=prc ) )

   DO k = 1, SIZE( beam%amp, 2 )
   	DO j = 1, SIZE( beam%amp, 1 )
   		beam%amp(j,k) = const * beam%amp(j,k)
   	END DO
   END DO

   END SUBROUTINE ApplyPiston



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE PistonTiltFocus( beam, piston, xTilt, yTilt, focus, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Explicitly applies piston, tilt and focus to the array of field values
   ! Piston input is OPD in meters, tilts are angles in radians, and focus is
   ! 1/radius of curvature. Tilt and focus aberrations are zero at x = xc,
   ! y = yc (default values of xc and yc are zero).
   ! origin.
   !     19 Oct  '00
   !     16 July '01    Revised to call QuadraticPhase
   !     12 Dec  '01    Revised to not require x() and y() work arrays
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines

   IMPLICIT NONE

   TYPE(wavefront),    INTENT(IN OUT) :: beam
   REAL(pr),           INTENT(IN)     :: piston
   REAL(pr),           INTENT(IN)     :: xTilt, yTilt
   REAL(pr),           INTENT(IN)     :: focus
   REAL(pr), OPTIONAL, INTENT(IN)     :: xc, yc

   INTRINSIC PRESENT

   ! Local variabales
   REAL(pr) :: pC, xT, yT, fC
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   pC = TWOPI*piston/beam%wavelength
   xT = TWOPI*xTilt/beam%wavelength
   yT = TWOPI*yTilt/beam%wavelength
   fC = PI*focus/beam%wavelength

   IF ( PRESENT( xc ) ) THEN
   	pC = pC - xc*( xT - (xc*fC) )
      xT = xT - 2.0_pr*(xc*fC)
   END IF
   IF ( PRESENT( yc ) ) THEN
   	pC = pC - yc*( yT - (yc*fC) )
      yT = yT - 2.0_pr*(yc*fC)
   END IF

   CALL QuadraticPhase( beam%amp, beam%x, beam%y, pC, xT, yT, fC )

   END SUBROUTINE PistonTiltFocus



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ApplyOPD( OPD, beam  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !     25 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN)     :: OPD
   TYPE(wavefront),          INTENT(IN OUT) :: beam

   INTRINSIC SIZE, EXP

   ! Local variables
   INTEGER  :: Ny, Nx
   INTEGER  :: j,  k
   REAL(pr) :: const
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( beam%amp, 2 )
   Ny = SIZE( beam%amp, 1 )

   IF ( SIZE( OPD, 1 ) /= Ny  .OR.  SIZE( OPD, 2 ) /= Nx ) &
      CALL ErrorExit( "ApplyOPD: array dimensions do not match" )

   const = TWOPI/beam%wavelength
   DO k = 1,Nx
   	DO j = 1,Ny
      	beam%amp(j,k) = beam%amp(j,k) * EXP( i*const*OPD(j,k) )
  		END DO
   END DO

   END SUBROUTINE ApplyOPD



	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	 SUBROUTINE ClearZernikeAberration( zAb )
	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	 ! Set values of zernikeAberration object to zero

	IMPLICIT NONE

   TYPE(zernikeAberration), INTENT(IN OUT) :: zAb

	! Local variables (none)
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	zAb%refDiam    = 0.0_pr
	zAb%xo         = 0.0_pr
	zAb%yo         = 0.0_pr
	zAb%modeNum    = 0
	zAb%coeffValue = 0.0_pr
	zAb%isRMS      = .FALSE.
	zAb%maxModes   = 0

	END SUBROUTINE ClearZernikeAberration



	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	 SUBROUTINE SetMultipleZernikeAberr( zAb, modes, coeffs, diam, xo, yo, RMS )
	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	 ! Set values of zernikeAberration object
	 !
	 ! diam		Reference diameter for Zernike polynomials
	 ! modes()	Array of mode numbers
	 ! coeffs()	Array of coefficient values (OPDs in SI units)
	 ! RMS		OPTIONAL logical flag for coefficients being RMS values,
	 !				rather than peak-to-valley (default is P-V)

	IMPLICIT NONE

	TYPE( zernikeAberration ), INTENT(IN OUT) :: zAb

	INTEGER,  DIMENSION(:),    INTENT(IN)     :: modes
	REAL(pr), DIMENSION(:),    INTENT(IN)     :: coeffs
	REAL(pr),                  INTENT(IN)     :: diam
	REAL(pr), OPTIONAL,        INTENT(IN)     :: xo, yo
	LOGICAL,  OPTIONAL,        INTENT(IN)     :: RMS

	INTRINSIC SIZE, PRESENT

	! Local variables:
	INTEGER :: modeIndex
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	IF ( SIZE( modes ) /= SIZE( coeffs )  .OR.  &
		  SIZE( modes ) >  MAX_ZERNIKE_MODES     ) THEN
		CALL ErrorExit( "SetMultipleZernikeAberr: array size error" )
	END IF

	IF ( diam <= 0.0_pr )  &
		CALL ErrorExit( "SetMultipleZernikeAberr: unphysical reference diameter" )

	CALL ClearZernikeAberration( zAb )

	DO modeIndex = SIZE( modes ), 1, -1
		IF ( modes(modeIndex) > 0 ) THEN
			zAb%maxModes = modeIndex
			EXIT
		END IF
	END DO

	zAb%refDiam                    = diam
	zAb%modeNum(1:zAb%maxModes)    = modes(1:zAb%maxModes)
	zAb%coeffValue(1:zAb%maxModes) = coeffs(1:zAb%maxModes)

	IF ( PRESENT( xo ) ) zAb%xo = xo
	IF ( PRESENT( yo ) ) zAb%yo = yo

	IF ( PRESENT( RMS ) ) zAb%isRMS = RMS

	END SUBROUTINE SetMultipleZernikeAberr


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	SUBROUTINE SetOneZernikeAberr( zAb, mode, coeff, diam, xo, yo, RMS )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	! Set values of zernikeAberration object
	!
	! diam		Reference diameter for Zernike polynomials
	! mode		Mode number
	! coeff		Coefficient value (OPD in SI units)
	! RMS		OPTIONAL logical flag for coefficient being RMS value,
	!				rather than peak-to-valley (default is P-V)

	IMPLICIT NONE

	TYPE( zernikeAberration ), INTENT(IN OUT) :: zAb

	INTEGER,  INTENT(IN)           :: mode
	REAL(pr), INTENT(IN)           :: coeff
	REAL(pr), INTENT(IN)           :: diam
	REAL(pr), INTENT(IN), OPTIONAL :: xo, yo
	LOGICAL,  INTENT(IN), OPTIONAL :: RMS

	INTRINSIC SIZE, PRESENT

	! Local variables:
	INTEGER :: modeIndex
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	IF ( diam <= 0.0_pr )  &
		CALL ErrorExit( "SetOneZernikeAberr: unphysical reference diameter" )

	CALL ClearZernikeAberration( zAb )

	zAb%maxModes      = 1
	zAb%refDiam       = diam
	zAb%modeNum(1)    = mode
	zAb%coeffValue(1) = coeff

	IF ( PRESENT( xo ) ) zAb%xo = xo
	IF ( PRESENT( yo ) ) zAb%yo = yo

	IF ( PRESENT( RMS ) ) zAb%isRMS = RMS

	END SUBROUTINE SetOneZernikeAberr



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ApplyOneZernike( beam, refDiam, mode, OPDcoeff, RMS, xo, yo  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Applies a Zernike phase to the input wavefront. Mode designations adopted
   ! from "Engineering & Laboratory Notes" from the August 1994 issue of
   ! Optics & Photonics News, by V.N. Mahajan.
   ! See also: R.J. Noll, JOSA vol.66, 207-211 (1976)
   !
   !  refDiam  Diameter of the circle within which the Zernikes are defined
   !  mode  	An index in the range 1 to 36 that specifies the Zernike mode
   !  OPDcoeff Magnitude of the mode, in meters of OPD
   !  RMS      An optional input; .TRUE. if OPDcoeff is an RMS value;
   !           otherwise, the coefficient is taken to be peak-to-valley
   !
   !  Note: although the Zernike polynomials are orthogonal only within the
   !        unit circle, this routine will calculate and apply phase values
   !        over the entire x-y plane. The value of refDiam must be chosen
   !        appropriately by the user.
   !
   !     17 Sept '01    Renamed ApplyOneZernike
   !     18 Dec  '01    Calls SUBROUTINE Zernike; uses ApplyOPD
   !      8 Oct  '02    Initializes OPD array to zero; OPTIONAL xo, yo
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Math_Routines, ONLY : Zernike

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr),            INTENT(IN) :: refDiam
   INTEGER,             INTENT(IN) :: mode
   REAL(pr),            INTENT(IN) :: OPDcoeff
   LOGICAL,  OPTIONAL,  INTENT(IN) :: RMS
   REAL(pr), OPTIONAL,  INTENT(IN) :: xo, yo

   INTRINSIC SIZE, PRESENT

   ! Local variables
   LOGICAL  :: isRMS
   REAL(pr), DIMENSION(SIZE( beam%y ),SIZE( beam%x )) :: OPD
   REAL(pr),                DIMENSION(SIZE( beam%x )) :: x
   REAL(pr),                DIMENSION(SIZE( beam%y )) :: y
   REAL(pr) :: xc, yc
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( mode <= 0 ) RETURN

   IF ( refDiam <= 0.0_pr ) CALL ErrorExit( "ApplyOneZernike: invalid refDiam" )

   xc = 0.0_pr
   yc = 0.0_pr
   IF ( PRESENT( xo ) ) xc = xo
   IF ( PRESENT( yo ) ) yc = yo

   x = (beam%x - xo)/(refDiam/2.0_pr)
   y = (beam%y - yo)/(refDiam/2.0_pr)

   IF ( PRESENT( RMS ) ) THEN
      isRMS = RMS
   ELSE
      isRMS = .FALSE.
   END IF

   OPD = 0.0_pr

   CALL Zernike( OPD, x, y, mode, OPDcoeff, isRMS )
   CALL ApplyOPD( OPD, beam )

   END SUBROUTINE ApplyOneZernike



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ApplyZernikeSum( beam, refDiam, modes, OPDcoeffs, RMS, xo, yo  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Applies a phase which is a sum of Zernikes.
   !
   !  refDiam		Diameter of the circle within which the Zernikes are defined
   !  modes			An index in the range 1 to 36 that specifies the Zernike mode
   !  OPDcoeffs	Magnitude of the mode, in meters of OPD
   !  RMS			An optional input; .TRUE. if OPDcoeff is an RMS value;
   !					otherwise, the coefficient is taken to be peak-to-valley
   !
   !  Note: although the Zernike polynomials are orthogonal only within the
   !        unit circle, this routine will calculate and apply phase values
   !        over the entire x-y plane. The value of refDiam must be chosen
   !        appropriately by the user.
   !
   !     17 Sept '01
   !     18 Dec  '01    Calls SUBROUTINE Zernike; uses ApplyOPD
   !      8 Oct  '02    Initializes OPD array to zero; assumes Zernike routine
   !                    adds Zernike mode to initial array; Zernike centered
   !                    at (xo,yo).
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Math_Routines, ONLY : Zernike

   IMPLICIT NONE

   TYPE(wavefront),    INTENT(IN OUT) :: beam
   REAL(pr),               INTENT(IN) :: refDiam
   INTEGER,  DIMENSION(:), INTENT(IN) :: modes
   REAL(pr), DIMENSION(:), INTENT(IN) :: OPDcoeffs
   LOGICAL,  OPTIONAL,     INTENT(IN) :: RMS
   REAL(pr), OPTIONAL,     INTENT(IN) :: xo, yo

   INTRINSIC SIZE, PRESENT, ABS, EPSILON

   ! Local variables
   LOGICAL  :: isRMS
   INTEGER  :: modeIndex
   REAL(pr), DIMENSION(SIZE( beam%y ),SIZE( beam%x )) :: OPD
   REAL(pr),                DIMENSION(SIZE( beam%x )) :: x
   REAL(pr),                DIMENSION(SIZE( beam%y )) :: y
   REAL(pr) :: xc, yc
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( SIZE( OPDcoeffs ) /= SIZE( modes ) )  &
   	CALL ErrorExit( "ApplyZernikeSum: input parameter error" )

   IF ( refDiam <= 0.0_pr ) CALL ErrorExit( "ApplyZernikeSum: invalid refDiam" )

   xc = 0.0_pr
   yc = 0.0_pr
   IF ( PRESENT( xo ) ) xc = xo
   IF ( PRESENT( yo ) ) yc = yo

   x = (beam%x - xo)/(refDiam/2.0_pr)
   y = (beam%y - yo)/(refDiam/2.0_pr)

   IF ( PRESENT( RMS ) ) THEN
      isRMS = RMS
   ELSE
      isRMS = .FALSE.
   END IF

   OPD = 0.0_pr

   DO modeIndex = 1, SIZE( modes )
   	CALL Zernike( OPD, x, y, modes(modeIndex), OPDcoeffs(modeIndex), isRMS )
   END DO

   CALL ApplyOPD( OPD, beam )

   END SUBROUTINE ApplyZernikeSum



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ApplyAberration( beam, zAberr  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Applies a phase which is a sum of Zernikes.
   !      8 Oct  '02
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	IMPLICIT NONE

   TYPE(wavefront),         INTENT(IN OUT) :: beam
   TYPE(zernikeAberration), INTENT(IN)     :: zAberr

   ! Local variables
   INTEGER  :: maxModes
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   maxModes = zAberr%maxModes
	IF ( maxModes > 0 ) THEN
		CALL ApplyZernike( beam, zAberr%refDiam,                 &
      								 zAberr%modeNum(1:maxModes),     &
                               zAberr%coeffValue(1:maxModes),  &
                               zAberr%isRMS,                   &
                               zAberr%xo,                      &
                               zAberr%yo                       )
	END IF

   END SUBROUTINE ApplyAberration



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Propagate( beam, z, dxo, dyo, applyCurv )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Propagates wavefront. Uses an "un-normalized" FFT (no factor of 1/N for
   ! either forward or inverse transform).

   !  z           Propagation distance
   !  dxo, dyo    Sample spacing in final plane
   !  applyCurv   Optional flag for applying wavefront curvature to complex
   !					amplitude at final plane. Default is to NOT apply curvature.

   !     14 Nov  '00
   !      2 July '01    Optional input: applyCurv
   !      2 Aug  '01    Removed beam%curvature from call to ApplyCurvature
   !     30 Nov  '01    Use DO loops for scaling amplitude array (save storage)
   !     19 June '02    Corrected Fresnel calculation
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines, ONLY : CoordinateVector
   USE Discrete_Transforms

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT)       :: beam
   REAL(pr),        INTENT(IN)           :: z
   REAL(pr),        INTENT(IN), OPTIONAL :: dxo, dyo
   LOGICAL,         INTENT(IN), OPTIONAL :: applyCurv

   INTRINSIC ABS, SIZE, EPSILON, PRESENT, SIGN

   ! Local variables
   REAL(pr) :: scaleFactor
   INTEGER  :: Nx, Ny
   INTEGER  ::  k, j
   REAL(pr) :: Gx, Gy
   REAL(pr) :: lambda_z
   REAL(pr) :: Fren_x,  Fren_y
   REAL(pr) :: dxout,   dyout
   REAL(pr) :: lambda_z_eff
   REAL(pr) :: geomCurv
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( z ) <= EPSILON(0.0_pr) ) RETURN

   scaleFactor = 1.0_pr + ( z * beam%curvature )

   Nx = SIZE( beam%amp, 2 )
   Ny = SIZE( beam%amp, 1 )

   Gx = Nx * beam%dx         ! Grid widths in each direction
   Gy = Ny * beam%dy

   ! Calculate Fresnel numbers based on grid size and
   ! divergence or convergence of initial wavefront

   lambda_z = z * beam%wavelength

   Fren_x = ABS( scaleFactor*Gx**2/lambda_z )
   Fren_y = ABS( scaleFactor*Gy**2/lambda_z )

   IF ( Fren_x <= Nx  .AND.  Fren_y <= Ny ) THEN

      ! Small Fresnel number case ("far-field")

      dxout = ABS( lambda_z/Gx )		! Sample spacing if using an FFT (default)
      dyout = ABS( lambda_z/Gy )
      IF ( PRESENT( dxo )  .AND.  PRESENT( dyo ) ) THEN
         IF ( dxo <= dxout  .AND.  dyo <= dyout ) THEN
				dxout = dxo
				dyout = dyo
         ELSE
            CALL ErrorExit( "Propagate: far-field grid too large" )
         END IF
      END IF

      CALL ApplyCurvature( beam, scaleFactor/z )
      CALL Fraunhofer
      beam%curvature = 1.0_pr/z

   ELSE IF ( Fren_x > Nx  .AND.  Fren_y > Ny ) THEN

      ! Large Fresnel number case ("near-field")

      IF ( PRESENT( dxo )  .AND.  PRESENT( dyo ) ) THEN
         ! Grid scale factor must be the same for both transverse coordinates
         IF ( ABS( (dxo/beam%dx) - (dyo/beam%dy) ) <= GRID_CRIT ) THEN
            scaleFactor = SIGN( dxo/beam%dx, scaleFactor )
            geomCurv    = ( scaleFactor - 1.0_pr )/z
		      ! Put residual curvature onto complex field
				CALL ApplyCurvature( beam, beam%curvature - geomCurv )
				dxout = dxo
				dyout = dyo
         ELSE
            CALL ErrorExit( "Propagate: grid scale factor error" )
         END IF
      ELSE
         ! Grid scale factor is detemined by wavefront curvature
         geomCurv = beam%curvature
         dxout    = ABS( scaleFactor ) * beam%dx
         dyout    = ABS( scaleFactor ) * beam%dy
      END IF

      lambda_z_eff = lambda_z/scaleFactor

      CALL Fresnel
      beam%curvature = geomCurv/scaleFactor

   ELSE

      CALL ErrorExit( "Propagate: Fresnel number conflict" )

   END IF

   beam%dx = dxout
   beam%dy = dyout

   beam%x = CoordinateVector( Nx, beam%dx )
   beam%y = CoordinateVector( Ny, beam%dy )

   IF ( PRESENT( applyCurv ) ) THEN
      IF ( applyCurv ) CALL ApplyCurvature( beam )
   END IF

   CONTAINS

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE Fraunhofer
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      !  10 Nov  '00
      !  18 Dec  '01    Calls revised FFTshift
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      INTRINSIC ABS, SIZE

      ! Local variables
      COMPLEX(pr) :: const
      REAL(pr)    :: dfx, dfy
      INTEGER     :: j,   k
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      const = (-i/lambda_z) * (beam%dx*beam%dy)
      dfx   = dxout/ABS( lambda_z )
      dfy   = dyout/ABS( lambda_z )

      ! Use an FFT if possible; otherwise use Chirp-z

      IF ( ABS( Gx*dfx - 1.0_pr ) <= GRID_CRIT  .AND. &
           ABS( Gy*dfy - 1.0_pr ) <= GRID_CRIT )  THEN

         PRINT*, " Fraunhofer using FFT2D! "
         CALL FFTshift( beam%amp, inverse=.TRUE. )
         !CALL FFT2D( beam%amp )
         call DFFTW_EXECUTE(beam%fftw_plan2d_for)

         CALL FFTshift( beam%amp )     ! Zero spatial frequency at grid center
      ELSE

         PRINT*, " Fraunhofer using CZT2D! "
         CALL CZT2D( beam%amp, beam%dy, beam%dx, dfy, dfx, beam%cztplan )
      END IF

      DO k = 1, Nx
      	DO j = 1, Ny
      		beam%amp(j,k) = const * beam%amp(j,k)
      	END DO
      END DO

      END SUBROUTINE Fraunhofer

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE Fresnel
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      !  10 Nov  '00
      !  18 Dec  '01    Calls revised FFTshift
      !  19 June '02    Does not call CZT2D for diverging/converging coordinates
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      INTRINSIC SIZE, EXP, ABS

      ! Local variables
      LOGICAL      :: inv
      REAL(pr)     :: const
      REAL(pr)     :: dfx, dfy
      INTEGER      :: j, k
      COMPLEX(prc) :: propConst
      REAL(pr),     DIMENSION(SIZE(beam%x)) :: fx
      REAL(pr),     DIMENSION(SIZE(beam%y)) :: fy
      COMPLEX(prc), DIMENSION(SIZE(beam%x)) :: efx
      COMPLEX(prc), DIMENSION(SIZE(beam%y)) :: efy
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      print*, "FRESNEL DIFFRACTION!"

      dfx = 1/Gx
      dfy = 1/Gy
      fx = CoordinateVector( Nx, dfx )
      fy = CoordinateVector( Ny, dfy )

      ! Use FFT to calculate plane-wave spectrum of input wavefront

      CALL FFTshift( beam%amp, inverse=.TRUE. )
      !CALL FFT2D( beam%amp )
      call DFFTW_EXECUTE(beam%fftw_plan2d_for)
      CALL FFTshift( beam%amp )     		! Zero frequency at grid center

      propConst = -i*PI*lambda_z_eff

      efx = (beam%dx*dfx/scaleFactor)*EXP( propConst * fx**2 )
      efy =             (beam%dy*dfy)*EXP( propConst * fy**2 )
      DO k = 1,Nx
      	DO j = 1,Ny
         	beam%amp(j,k) = efy(j)*efx(k)*beam%amp(j,k)
      	END DO
      END DO

      IF ( scaleFactor < 0.0_pr ) THEN    ! Propagating through focus inverts
         inv = .FALSE.                    ! the wavefront
      ELSE
         inv = .TRUE.
      END IF

      CALL FFTshift( beam%amp, inverse=.TRUE. )
      !CALL FFT2D( beam%amp, inverse = inv )
      call DFFTW_EXECUTE(beam%fftw_plan2d_rev)

      CALL FFTshift( beam%amp )

      END SUBROUTINE Fresnel

   END SUBROUTINE Propagate



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Apodize( beam, apodizer, outerDiam, innerDiam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Applies a circularly symmetric apodization, the apodizer type and outer
   ! diameter must always be specified; the inner diameter and the location of
   ! the apodizer center (xc,yc) have default values of zero in the apodization
   ! routines. The outer diameter is the circle with center (xc,yc) on which
   ! the wavefront amplitude is reduced to 1/e of its initial value for Gaussian
   ! taper, and to zero for cosine or quadratic taper.

   !     11 Mar  '02
   !     20 Sept '02    Corrected comment only
   !     10 Oct  '02    Added Quadratic apodization
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Clipping_Routines

   IMPLICIT NONE

   TYPE(wavefront),  INTENT(IN OUT)       :: beam
   CHARACTER(LEN=*), INTENT(IN)           :: apodizer
   REAL(pr),         INTENT(IN)           :: outerDiam
   REAL(pr),         INTENT(IN), OPTIONAL :: innerDiam
   REAL(pr),         INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC PRESENT

   ! Local variables
   REAL(pr) :: innerD
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   innerD = 0.0_pr
   IF ( PRESENT( innerDiam ) ) innerD = innerDiam

   SELECT CASE ( apodizer )
      CASE ( "gaussian", "Gaussian", "GAUSSIAN" )
         CALL GaussianCircularApodization( beam%amp, beam%x, beam%y, &
                                           outerDiam, innerD, xc, yc )
      CASE ( "cosine", "Cosine", "COSINE" )
         CALL CosineCircularApodization( beam%amp, beam%x, beam%y, &
                                         outerDiam, innerD, xc, yc )
      CASE ( "quadratic", "Quadratic", "QUADRATIC" )
         CALL QuadraticCircularApodization( beam%amp, beam%x, beam%y, &
														  outerDiam, innerD, xc, yc )
      CASE DEFAULT
         CALL ErrorExit( "ApodizeWavefront: apodizer not correctly specified" )
   END SELECT

   END SUBROUTINE Apodize



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ClipCirc( beam, diam, xc, yc, obsc, tiltAngle, tiltOrientation )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calls for clipping -- circular aperture or obscuration. A tilted circular
   ! aperture or obscuration is rotated about an axis in the x-y plane which
   ! makes an angle "tiltOrientation" with respect to the x-axis; the amount
   ! of tilt is "tiltAngle."

   !     26 Oct  '00
   !      6 July '01    Tests for zero diam
   !     20 Sept '02    Aperture may be tilted
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Clipping_Routines

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT)       :: beam
   REAL(pr),        INTENT(IN)           :: diam
   REAL(pr),        INTENT(IN), OPTIONAL :: xc, yc
   LOGICAL,         INTENT(IN), OPTIONAL :: obsc
   REAL(pr),        INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

   INTRINSIC EPSILON, PRESENT, ABS, COS

   ! Local variables
   LOGICAL  :: isObscuration
   LOGICAL  :: isTilted
   REAL(pr) :: minorDiam
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( diam < EPSILON( diam ) ) RETURN

   isObscuration = .FALSE.
   IF ( PRESENT( obsc ) ) isObscuration = obsc

   isTilted = .FALSE.
   IF ( PRESENT( tiltAngle )  .AND.  PRESENT( tiltOrientation ) ) THEN
      IF ( ABS( tiltAngle ) > EPSILON( tiltAngle ) ) THEN
      	isTilted = .TRUE.
      END IF
   END IF

   IF ( isTilted ) THEN

      minorDiam = diam * COS( tiltAngle )

      IF ( isObscuration ) THEN
         CALL EllipticalObsc( beam%amp, beam%x, beam%y, diam, minorDiam,  &
      							   tiltOrientation, xc, yc )
      ELSE
         CALL EllipticalClip( beam%amp, beam%x, beam%y, diam, minorDiam,  &
      							   tiltOrientation, xc, yc )
      END IF

   ELSE

      IF ( isObscuration ) THEN
         CALL CircularObsc( beam%amp, beam%x, beam%y, diam, xc, yc )
      ELSE
         CALL CircularClip( beam%amp, beam%x, beam%y, diam, xc, yc )
      END IF

   END IF

   END SUBROUTINE ClipCirc



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ClipRect( beam, length, width, xc, yc, angle, obsc,  &
   							tiltAngle, tiltOrientation )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Rectangular aperture or obscuration. Angle is the counterclockwise rotation
   ! of the rectangle about the optical axis; zero angle leaves "length" along
   ! x-axis, "width" along y. If the aperture is tilted, the angle between the
   ! x-axis and the rotation axis for the tilt is specified by tiltOrientation,
   ! while tiltAngle is the magnitude of the tilt.

   !     26 Oct  '00
   !      6 July '01    Tests for zero-size rectangle
   !     19 July '01    Clears "rectangle" after use
   !     20 Sept '02    Handles tilted aperture/obscuration
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons
   USE Clipping_Routines

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT)       :: beam
   REAL(pr),        INTENT(IN)           :: length, width
   REAL(pr),        INTENT(IN), OPTIONAL :: xc, yc
   REAL(pr),        INTENT(IN), OPTIONAL :: angle
   LOGICAL,         INTENT(IN), OPTIONAL :: obsc
   REAL(pr),        INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

   INTRINSIC EPSILON, PRESENT, ABS, COS, SIN

   ! Local variables
   TYPE(polygon) :: rectangle
   LOGICAL       :: isObscuration
   REAL(pr), DIMENSION(3) :: axisVector
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF (  length < EPSILON( length )  .OR.  width < EPSILON( width ) ) RETURN

   isObscuration = .FALSE.
   IF ( PRESENT( obsc ) ) isObscuration = obsc

   rectangle = MakeRectangle( length, width, angle, xc, yc )

   IF ( PRESENT( tiltAngle )  .AND.  PRESENT( tiltOrientation ) ) THEN
      IF ( ABS( tiltAngle ) > EPSILON( tiltAngle ) ) THEN
         axisVector = (/ COS(tiltOrientation), SIN(tiltOrientation), 0.0_pr /)
         CALL RotatePolygon_3d( rectangle, tiltAngle, axisVector )
      END IF
   END IF

   IF ( isObscuration ) THEN
      CALL PolygonObsc( beam%amp, beam%x, beam%y, rectangle )
   ELSE
      CALL PolygonClip( beam%amp, beam%x, beam%y, rectangle )
   END IF

   CALL ClearPolygon( rectangle )

   END SUBROUTINE ClipRect



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ClipPoly( beam, poly, xc, yc, angle, obsc,  &
   							tiltAngle, tiltOrientation )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clipping by polygonal aperture or obscuration. xc and yc are shifts of the
   ! input polygon, and angle is a counterclockwise rotation. The input polygon
   ! is rotated about (x=0,y=0) before being shifted. If the aperture is tilted,
   ! the angle between the x-axis and the rotation axis for the tilt is specified
   ! by tiltOrientation, while tiltAngle is the magnitude of the tilt.

   !     26 Oct  '00
   !      6 July '01    Tests for zero polygon area
   !     24 Sept '01    Rotation implemented
   !     20 Sept '02    Aperture tilt implemented
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons
   USE Clipping_Routines

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT)       :: beam
   TYPE(polygon),   INTENT(IN)           :: poly
   REAL(pr),        INTENT(IN), OPTIONAL :: xc, yc
   REAL(pr),        INTENT(IN), OPTIONAL :: angle
   LOGICAL,         INTENT(IN), OPTIONAL :: obsc
   REAL(pr),        INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

   INTRINSIC EPSILON, PRESENT, ABS, COS, SIN

   ! Local variables
   TYPE(polygon) :: shiftedPoly
   LOGICAL       :: isObscuration
   REAL(pr)      :: xShift, yShift
   REAL(pr), DIMENSION(3) :: axisVector
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( PolygonArea( poly ) < EPSILON( 0.0_pr ) ) RETURN

   shiftedPoly = poly

   IF ( PRESENT( angle ) ) CALL RotatePolygon( shiftedPoly, angle )

   xShift = 0.0_pr
   yShift = 0.0_pr
   IF ( PRESENT( xc ) ) xShift = xc
   IF ( PRESENT( yc ) ) yShift = yc
   CALL ShiftPolygon( shiftedPoly, xShift, yShift )

   IF ( PRESENT( tiltAngle )  .AND.  PRESENT( tiltOrientation ) ) THEN
      IF ( ABS( tiltAngle ) > EPSILON( tiltAngle) ) THEN
         axisVector = (/ COS(tiltOrientation), SIN(tiltOrientation), 0.0_pr /)
         CALL RotatePolygon_3d( shiftedPoly, tiltAngle, axisVector )
      END IF
   END IF

   isObscuration = .FALSE.
   IF ( PRESENT( obsc ) ) isObscuration = obsc

   IF ( isObscuration ) THEN
      CALL PolygonObsc( beam%amp, beam%x, beam%y, shiftedPoly )
   ELSE
      CALL PolygonClip( beam%amp, beam%x, beam%y, shiftedPoly )
   END IF

   CALL ClearPolygon( shiftedPoly )

   END SUBROUTINE ClipPoly



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SpiderObsc( beam, numVanes, width, startAngle, xc, yc,  &
   							  tiltAngle, tiltOrientation )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Obscuration due to spider support. The vanes will be equally spaced in
   ! angle, radiating out from the point (xc,yc). If the spiders are tilted,
   ! the angle between the x-axis and the rotation axis for the tilt is
   ! specified by tiltOrientation, while tiltAngle is the magnitude of the tilt.
   !
   !  numVanes   The number of vanes (usually 3 or 4; at least 1)
   !  width		   Width of the vanes (narrow dimension)
   !  startAngle  Angle between x-axis and the first vane; default is zero
   !  xc, yc      Coordinates of point where vanes meet; default is (0,0)

   !     24 Sept '01
   !      1 Oct  '01    Corrected arguments for MakeRectangle
   !     20 Sept '02    Tilt implemented
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons
   USE Clipping_Routines
   USE Constants

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT)       :: beam
   INTEGER,         INTENT(IN)           :: numVanes
   REAL(pr),        INTENT(IN)           :: width
   REAL(pr),        INTENT(IN), OPTIONAL :: startAngle
   REAL(pr),        INTENT(IN), OPTIONAL :: xc, yc
   REAL(pr),        INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

   INTRINSIC SIZE, PRESENT, SQRT, MAX, COS, SIN

   ! Local variables
   TYPE(polygon) :: vane, tiltedVane
   INTEGER       :: n
   REAL(pr)      :: length
   INTEGER       :: Nx, Ny
   REAL(pr)      :: xCenter, yCenter
   REAL(pr)      :: angle, theta, rotAngle
   LOGICAL       :: isTilted
   REAL(pr), DIMENSION(3) :: axisVector
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( beam%x )
   Ny = SIZE( beam%y )

   angle = 0.0_pr
   IF ( PRESENT( startAngle ) ) angle = startAngle

   isTilted = .FALSE.
   IF ( PRESENT( tiltAngle )  .AND.  PRESENT( tiltOrientation ) ) THEN
      IF ( ABS( tiltAngle ) > EPSILON( tiltAngle ) ) THEN
         isTilted   = .TRUE.
         axisVector = (/ COS(tiltOrientation), SIN(tiltOrientation), 0.0_pr /)
      END IF
   END IF

   ! Angular spacing of vanes
   theta = TWOPI/REAL( numVanes, KIND=pr )

   xCenter = 0.0_pr
   yCenter = 0.0_pr
   IF ( PRESENT( xc ) ) xCenter = xc
   IF ( PRESENT( yc ) ) yCenter = yc

   ! Make a vane long enough to reach to edge of beam grid
   length = SQRT( ( MAX( (beam%x(1)-xCenter)**2, (beam%x(Nx)-xCenter)**2 )  + &
                    MAX( (beam%y(1)-yCenter)**2, (beam%y(Ny)-yCenter)**2 ) )  )

   ! Create a vane with it's end at (xc,yc) and parallel to the x-axis.
   vane = MakeRectangle( length, width, 0.0_pr, 0.5_pr*length+xCenter, yCenter )

   ! Rotate the vane about (xc,yc) and clip the beam; repeat for each vane
   rotAngle = angle
   DO n = 1,numVanes
      CALL RotatePolygon( vane, rotAngle, xCenter, yCenter )
      IF ( isTilted ) THEN
         tiltedVane = vane
      	CALL RotatePolygon_3d( tiltedVane, tiltAngle, axisVector )
      	CALL PolygonObsc( beam%amp, beam%x, beam%y, tiltedVane )
      ELSE
      	CALL PolygonObsc( beam%amp, beam%x, beam%y, vane )
      END IF
      rotAngle = theta
   END DO

   CALL ClearPolygon( vane )
   CALL ClearPolygon( tiltedVane )

   END SUBROUTINE SpiderObsc



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE InvertWavefrontX( beam  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! The input wavefront is right-left inverted (as when reflected from mirror)
   !      6 Feb  '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : ReflectAboutYaxis

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam

   ! Local variables: (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   CALL ReflectAboutYaxis( beam%amp )

   END SUBROUTINE InvertWavefrontX



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE InvertWavefrontY( beam  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! The input wavefront is up-down inverted (as when reflected from mirror)
   !      6 Feb  '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : ReflectAboutXaxis

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam

   ! Local variables: (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   CALL ReflectAboutXaxis( beam%amp )

   END SUBROUTINE InvertWavefrontY



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE InvertWavefront( beam, xc, yc  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! The input wavefront is inverted through point (xc,yc); default is (0,0)

   !      6 Jan  '01
   !      2 Feb  '01    Calls revised ReflectThroughOrigin
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : ReflectThroughOrigin

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT)       :: beam
   REAL(pr),        INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC PRESENT

   ! Local variables: (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( PRESENT( xc )  .AND. PRESENT( yc ) ) THEN
      CALL ReflectThroughOrigin( beam%amp, xc/beam%dx, yc/beam%dy )
   ELSE
      CALL ReflectThroughOrigin( beam%amp )
   END IF

   END SUBROUTINE InvertWavefront



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ReflectFromCornerCube( cc, beam  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! The input wavefront is passed to CornerCubeReflect

   !     26 Oct  '00
   !      2 Feb  '01    Calls revised ReflectThroughOrigin
   !     19 July '01    Clears clearAp after use
   !      4 Dec  '01    Rewritten to call revised CornerCubeReflect
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Corner_Cubes

   IMPLICIT NONE

   TYPE(CornerCube), INTENT(IN)     :: cc
   TYPE(wavefront),  INTENT(IN OUT) :: beam

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   CALL CornerCubeReflect( cc, beam%amp, beam%dx, beam%x, beam%dy, beam%y,  &
   									 beam%wavelength )

   END SUBROUTINE ReflectFromCornerCube



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION GetWavefrontIntensity( beam, Nx, Ny ) RESULT( intensity )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns an array containing the intensity of the input beam. The Nx and Ny
   ! samples are centered on the center of the beam grid.
   !     19 Oct  '00
   !      7 June '01    Number of output samples (Nx,Ny) is required input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   INTEGER,         INTENT(IN) :: Nx, Ny

   REAL(pr),  DIMENSION(Ny,Nx) :: intensity

   INTRINSIC SIZE, PRESENT, REAL, AIMAG

   ! Local variables
   INTEGER :: j,     k
   INTEGER :: rows,  cols
   INTEGER :: jDiff, kDiff
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   rows = SIZE(beam%amp,1)
   cols = SIZE(beam%amp,2)

   IF ( Ny > rows  .OR.  Nx > cols ) THEN
      CALL ErrorExit( "Nx or Ny too large in GetWavefrontIntensity" )
   END IF

   jDiff = (rows/2) - (Ny/2)
   kDiff = (cols/2) - (Nx/2)

   DO k = 1,Nx
   	DO j = 1,Ny
      	intensity(j,k) = REAL( beam%amp(j+jDiff,k+kDiff) )**2 + &
                       	 AIMAG( beam%amp(j+jDiff,k+kDiff) )**2
   	END DO
   END DO

   END FUNCTION GetWavefrontIntensity



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION GetWavefrontAmp( beam, Nx, Ny ) RESULT( amplitude )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns an array containing the complex amplitude of the input beam.
   ! The Nx and Ny samples are centered on the center of the beam grid.
   !     20 Nov  '00
   !      7 June '01    Number of output samples (Nx,Ny) is required input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   INTEGER,         INTENT(IN) :: Nx, Ny

   COMPLEX(prc), DIMENSION(Ny,Nx) :: amplitude

   INTRINSIC SIZE

   ! Local variables
   INTEGER j,     k
   INTEGER rows,  cols
   INTEGER jDiff, kDiff
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   rows = SIZE(beam%amp,1)
   cols = SIZE(beam%amp,2)

   IF ( Ny > rows  .OR.  Nx > cols ) THEN
      CALL ErrorExit( "Nx or Ny too large in GetWavefrontAmp" )
   END IF

   jDiff = (rows/2) - (Ny/2)
   kDiff = (cols/2) - (Nx/2)

   DO k = 1,Nx
   	DO j = 1,Ny
         amplitude(j,k) = beam%amp(j+jDiff,k+kDiff)
      END DO
   END DO

   END FUNCTION GetWavefrontAmp



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION GetWavefrontOPD( beam, Nx, Ny ) RESULT( OPD )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns an array containing the OPD of the input beam. The Nx and Ny
   ! samples are centered on the center of the beam grid.
   !      7 June '01    Number of output samples (Nx,Ny) is required input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Math_Routines, ONLY : PhaseAngle

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   INTEGER,         INTENT(IN) :: Nx, Ny

   REAL(pr),  DIMENSION(Ny,Nx) :: OPD

   INTRINSIC SIZE

   ! Local variables
   INTEGER  :: j,     k
   INTEGER  :: rows,  cols
   INTEGER  :: jDiff, kDiff
   REAL(pr) :: const
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   rows = SIZE(beam%amp,1)
   cols = SIZE(beam%amp,2)

   IF ( Ny > rows  .OR.  Nx > cols ) THEN
      CALL ErrorExit( "Nx or Ny too large in GetWavefrontOPD" )
   END IF

   jDiff = (rows/2) - (Ny/2)
   kDiff = (cols/2) - (Nx/2)

   const = beam%wavelength/TWOPI

   DO k = 1,Nx
   	DO j = 1,Ny
         OPD(j,k) = const*PhaseAngle( beam%amp(j+jDiff,k+kDiff) )
      END DO
   END DO

   END FUNCTION GetWavefrontOPD



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GetWavefrontSection_XY( beam, x1, x2, y1, y2, section )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns a wavefront section that contains the portion of the
   ! input wavefront that lies in the region (x1 <= x <= x2, y1 <= y <= y2).

   !     25 Sept '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   REAL(pr),        INTENT(IN) :: x1, x2, y1, y2

   TYPE(wavefront), INTENT(IN OUT) :: section

   INTRINSIC SIZE

   ! Local variables
   INTEGER :: Ny,  Nx
   INTEGER ::  j,  k
   INTEGER :: Nys, Nxs
   INTEGER :: j1, j2, k1, k2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE(beam%amp,1)
   Nx = SIZE(beam%amp,2)

   ! Find largest j1 in beam coordinate array for which y(j1) <= y1,
   ! and the smallest j2 such that y(j2) >= y2
   j1 = 0
   j2 = 0
   DO j = Ny,1,-1
      IF ( beam%y(j) <= y1 ) THEN
         j1 = j
         EXIT
      END IF
   END DO
   DO j = 1,Ny
      IF ( beam%y(j) >= y2 ) THEN
         j2 = j
         EXIT
      END IF
   END DO
   IF ( j1 == 0  .OR.  j2 == 0 )  &
      CALL ErrorExit( "GetWavefrontSection_XY: can't find acceptable y-range" )

   ! Find corresponding x-interval
   k1 = 0
   k2 = 0
   DO k = Nx,1,-1
      IF ( beam%x(k) <= x1 ) THEN
         k1 = k
         EXIT
      END IF
   END DO
   DO k = 1,Nx
      IF ( beam%x(k) >= x2 ) THEN
         k2 = k
         EXIT
      END IF
   END DO
   IF ( k1 == 0  .OR.  k2 == 0 )  &
      CALL ErrorExit( "GetWavefrontSection_XY: can't find acceptable x-range" )

   ! Number of rows and columns for output wavefront section:
   Nys = j2-j1+1
   Nxs = k2-k1+1

   CALL TestWavefrontAllocation( section, Nxs, Nys )
   section%dx = beam%dx
   section%dy = beam%dy
   section%wavelength = beam%wavelength
   section%curvature  = beam%curvature
   section%x(:) = beam%x(k1:k2)
   section%y(:) = beam%y(j1:j2)
   section%amp(:,:) = beam%amp(j1:j2,k1:k2)

   END SUBROUTINE GetWavefrontSection_XY



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GetWavefrontSection_JK( beam, j1, j2, k1, k2, section )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns a wavefront section that contains the portion of the
   ! input wavefront with rows and columns (j1 <= j <= j2, k1 <= k <= k2),
   ! where the amplitude array rows correspond to y(j), and columns to x(k).

   !     25 Sept '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   INTEGER,         INTENT(IN) :: j1, j2, k1, k2

   TYPE(wavefront), INTENT(IN OUT) :: section

   INTRINSIC SIZE

   ! Local variables
   INTEGER :: Ny,  Nx
   INTEGER :: Nys, Nxs
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE(beam%amp,1)
   Nx = SIZE(beam%amp,2)

   IF ( j1 < 1  .OR.  j1 > Ny  .OR.  j2 < 1  .OR.  j2 > Ny )  &
      CALL ErrorExit( "GetWavefrontSection_JK: j1, j2 not consistent with beam" )

   IF ( k1 < 1  .OR.  k1 > Nx  .OR.  k2 < 1  .OR.  k2 > Nx )  &
      CALL ErrorExit( "GetWavefrontSection_JK: k1, k2 not consistent with beam" )

   ! Number of rows and columns for output wavefront section:
   Nys = j2-j1+1
   Nxs = k2-k1+1

   CALL TestWavefrontAllocation( section, Nxs, Nys )
   section%dx = beam%dx
   section%dy = beam%dy
   section%wavelength = beam%wavelength
   section%curvature  = beam%curvature
   section%x(:) = beam%x(k1:k2)
   section%y(:) = beam%y(j1:j2)
   section%amp(:,:) = beam%amp(j1:j2,k1:k2)

   END SUBROUTINE GetWavefrontSection_JK



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE AddWavefrontSection( beam1, beam2 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Adds the complex values of two input wavefronts; beam2 is a wavefront
   ! section (i.e., has same spacing as beam1, but fewer samples).
   ! If the input wavefronts have the same number of sample points, then this
   ! routine should produce the same output as WavefrontSum().
   !     25 Sept '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines, ONLY : QuadraticPhase

   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam1  ! "Principal" wavefront
   TYPE(wavefront), INTENT(IN)     :: beam2  ! Wavefront section

   INTRINSIC ABS, SIZE, MAX

   ! Local variables
   INTEGER  :: Ny1, Nx1, Ny2, Nx2
   REAL(pr) :: deltaLambda
   INTEGER  :: j, k
   INTEGER  :: j1, j2, k1, k2
   REAL(pr) :: curvCoeff, piston, tilt
   COMPLEX(prc), DIMENSION(SIZE(beam2%amp,1),SIZE(beam2%amp,2)) :: tempAmp
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


   deltaLambda = beam1%wavelength - beam2%wavelength

	IF ( ( ABS( beam1%dx - beam2%dx ) > beam1%dx*GRID_CRIT )  .OR.  &
		  ( ABS( beam1%dy - beam2%dy ) > beam1%dy*GRID_CRIT )  ) THEN

      CALL ErrorExit( "AddWavefrontSection: must have equal sample spacing" )

   ELSE IF ( ABS( deltaLambda ) > beam1%wavelength*WAVELENGTH_CRIT ) THEN

      CALL ErrorExit( "AddWavefrontSection: must have equal wavelengths" )

   ELSE

      Ny1 = SIZE( beam1%amp, 1 )
      Nx1 = SIZE( beam1%amp, 2 )
      Ny2 = SIZE( beam2%amp, 1 )
      Nx2 = SIZE( beam2%amp, 2 )

      IF ( Nx1 >= Nx2  .AND.  Ny1 >= Ny2 ) THEN

         ! Beam#2 must be a section (i.e., must have the smaller grid)

         ! Apply differential curvature to amplitude of beam section
         tempAmp(:,:) = beam2%amp(:,:)
         piston = 0.0_pr
         tilt   = 0.0_pr
         curvCoeff = PI*(beam2%curvature - beam1%curvature)/beam2%wavelength
         CALL QuadraticPhase( tempAmp, beam2%x, beam2%y,    &
         							piston, tilt, tilt, curvCoeff )

         ! Find where beam2 fits into the grid of beam1
         k1 = 0
         k2 = 0
         DO k = 1,Nx1
         	IF ( ABS( beam1%x(k)-beam2%x(1) )   < beam1%dx*GRID_CRIT ) k1 = k
         	IF ( ABS( beam1%x(k)-beam2%x(Nx2) ) < beam1%dx*GRID_CRIT ) k2 = k
         END DO
         IF ( k1 == 0  .OR.  k2 == 0 )  &
            CALL ErrorExit( "AddWavefrontSection: x-index error" )
         j1 = 0
         j2 = 0
         DO j = 1,Ny1
         	IF ( ABS( beam1%y(j)-beam2%y(1) )   < beam1%dy*GRID_CRIT ) j1 = j
         	IF ( ABS( beam1%y(j)-beam2%y(Ny2) ) < beam1%dy*GRID_CRIT ) j2 = j
         END DO
         IF ( j1 == 0  .OR.  j2 == 0 )  &
            CALL ErrorExit( "AddWavefrontSection: y-index error" )

         ! Coherently add wavefront amplitudes at corresponding sample points
         beam1%amp(j1:j2,k1:k2) = beam1%amp(j1:j2,k1:k2) + tempAmp(:,:)

      ELSE
         CALL ErrorExit( "AddWavefrontSection: beam2 not a section of beam1" )
      END IF

   END IF

   END SUBROUTINE AddWavefrontSection



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Heterodyne( sigBeam, refBeam, OPD, power,   &
                          detShape, xSize, ySize, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the heterodyne signal associated with the interference of two
   ! input beams. The detector shape may be "circle," "square," "rectangle"
   ! or "infinite" (whole beam grid). For an "infinite" detector, the size and
   ! center parameters are not used. ySize parameter required only for
   ! a rectangular detector.

   !     31 Mar '01     Initial version, based on AverageOPD function
   !      7 Jun '01     Calls PhaseAngle, instead of ATAN2
   !      9 July '01    Shape "infinite" added
   !     16 July '01    Applies curvature to "hetAmp"
   !      2 Aug  '01    Corrected value of "curv" for quadratic phase
   !     30 Nov  '01    Use DO loops to save temp storage
   !     28 Feb  '02    Calls Quad2d for integrating over detector
   !      5 Mar  '02    "power" is magnitude of phasor
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Quadrature_Routines
   USE Math_Routines, ONLY : PhaseAngle
   USE Array_Routines

   IMPLICIT NONE

   TYPE(wavefront),   INTENT(IN)  :: sigBeam, refBeam

   REAL(pr),          INTENT(OUT) :: OPD, power

   CHARACTER( LEN=* ), INTENT(IN) :: detShape
   REAL(pr),           INTENT(IN), OPTIONAL :: xSize
   REAL(pr),           INTENT(IN), OPTIONAL :: ySize
   REAL(pr),           INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SIZE, ABS, CONJG

   ! Local variables
   COMPLEX(prc), DIMENSION(SIZE(sigBeam%amp,1),SIZE(sigBeam%amp,2)) :: hetAmp
   COMPLEX(prc) :: phasor
   REAL(pr)     :: piston, xTilt, yTilt, curv
   INTEGER      :: Nx, k, Ny, j
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( sigBeam%amp, 2 )
   Ny = SIZE( sigBeam%amp, 1 )

   IF ( SIZE(refBeam%amp,1) /= Ny                             .OR.  &
        SIZE(refBeam%amp,2) /= Nx                             .OR.  &
        ABS( refBeam%dx - sigBeam%dx ) > GRID_CRIT*sigBeam%dx .OR.  &
        ABS( refBeam%dy - sigBeam%dy ) > GRID_CRIT*sigBeam%dy )  THEN
      CALL ErrorExit( "Heterodyne: beam grids do not match" )
   END IF

   DO k = 1, Nx
   	DO j = 1, Ny
   		hetAmp(j,k) = sigBeam%amp(j,k) * CONJG( refBeam%amp(j,k) )
   	END DO
   END DO

   piston = 0.0_pr
   xTilt  = 0.0_pr
   yTilt  = 0.0_pr
   curv   = PI*(sigBeam%curvature-refBeam%curvature)/sigBeam%wavelength
   CALL QuadraticPhase( hetAmp, sigBeam%x, sigBeam%y,  &
                        piston, xTilt, yTilt, curv )

   phasor = Quad2d( hetAmp, sigBeam%x, sigBeam%y, detShape,  &
                    xSize, ySize, xc, yc )

   power = ABS( phasor )
   OPD   = (sigBeam%wavelength/TWOPI) * PhaseAngle( phasor )

   END SUBROUTINE Heterodyne



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION AverageOPD( beam, detShape, xSize, ySize, xc, yc ) RESULT( OPD )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculate the average OPD for a finite detector, whose shape may be
   ! "circle" "square", "rectangle" or "infinite."

   !     31 Oct  '00
   !      8 Nov  '00    Uses Quad2D routines
   !      7 Jun  '01    Calls PhaseAngle, instead of ATAN2
   !      9 July '01    Added option of "infinite" (whole grid) detector
   !      1 Mar  '02    Calls new Quad2d function; rearranged argument list
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Quadrature_Routines
   USE Math_Routines, ONLY : PhaseAngle

   IMPLICIT NONE

   TYPE(wavefront),    INTENT(IN) :: beam
   CHARACTER( LEN=* ), INTENT(IN) :: detShape
   REAL(pr),           INTENT(IN), OPTIONAL :: xSize  ! Width or diameter
   REAL(pr),           INTENT(IN), OPTIONAL :: ySize  ! Height (rectangle)
   REAL(pr),           INTENT(IN), OPTIONAL :: xc, yc ! Detector center

   REAL(pr) :: OPD

   ! Local variables
   COMPLEX(prc) :: phasor
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   phasor = Quad2d( beam%amp, beam%x, beam%y, detShape,  &
                    xSize, ySize, xc, yc )
   OPD    = (beam%wavelength/TWOPI) * PhaseAngle( phasor )

   END FUNCTION AverageOPD



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SaveWavefront( beam, unitNo )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Writes unformatted wavefront data to SCRATCH disk for temporary storage.
   ! If unit is already open, the unit is closed and then opened so the new
   ! wavefront can be written on it.
   !     19 Mar '01     R.S. Benson
   !     30 Oct '01     Tests for unitNo being already open
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN) :: beam
   INTEGER,         INTENT(IN) :: unitNo

   INTRINSIC SIZE

   ! Local variables
   INTEGER  :: Nx, Ny
   INTEGER  ::  k,  j
   REAL(pr) :: dx, dy
   REAL(pr) :: wavelength, curvature
   LOGICAL  :: isOpen
   INTEGER  :: ioFlag
   INTEGER, PARAMETER :: ok = 0
   CHARACTER(LEN=3)   :: ioUnit
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( beam%amp, 2 )
   Ny = SIZE( beam%amp, 1 )

   INQUIRE( unitNo, OPENED = isOpen )
   IF ( isOpen ) CLOSE( unitNo )             ! Get rid of old data on this unit
   OPEN( unitNo, STATUS   = "SCRATCH",     &
                 FORM     = "UNFORMATTED", &
                 ACTION   = "READWRITE",   &
                 POSITION = "REWIND",      &
                 IOSTAT   = ioFlag )
   IF ( ioFlag == ok ) THEN

      dx = beam%dx
      dy = beam%dy
      wavelength = beam%wavelength
      curvature  = beam%curvature
      WRITE( unitNo ) Nx, Ny
      WRITE( unitNo ) dx, dy, wavelength, curvature
      WRITE( unitNo ) ( ( beam%amp(j,k), j=1,Ny ), k=1,Nx )
      WRITE( unitNo ) ( beam%x(k), k=1,Nx ),  ( beam%y(j), j=1,Ny )
      REWIND( unitNo )

   ELSE
      WRITE( ioUnit, "(I3)" ) unitNo
      CALL ErrorExit( "SaveWavefront: error opening unit no. "//ioUnit )
   END IF

   END SUBROUTINE SaveWavefront



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE RetrieveWavefront( beam, unitNo )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Reads unformatted wavefront data from SCRATCH file.
   ! UNIT = unitNo should already have been opened.
   !     19 Mar '01     R.S. Benson
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   INTEGER,         INTENT(IN)     :: unitNo

   ! Local variables
   INTEGER  :: Nx, Ny
   INTEGER  ::  k,  j
   REAL(pr) :: dx, dy
   REAL(pr) :: wavelength, curvature
   INTEGER, PARAMETER :: ok = 0
   LOGICAL            :: isOpen
   CHARACTER(LEN=3)   :: ioUnit
   INTEGER            :: ioStatus
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   INQUIRE( unitNo, OPENED = isOpen )
   IF ( .NOT. isOpen ) THEN
      WRITE( ioUnit, "(I3)" ) unitNo
      CALL ErrorExit( "RetrieveWavefront: scratch unit "//ioUnit//" is not open" )
   END IF

   READ( unitNo, IOSTAT = ioStatus ) Nx, Ny

   IF ( ioStatus == ok ) THEN

      CALL TestWavefrontAllocation( beam, Nx, Ny )

      READ( unitNo ) dx, dy, wavelength, curvature
      READ( unitNo ) ( ( beam%amp(j,k), j=1,Ny ), k=1,Nx )
      READ( unitNo ) ( beam%x(k), k=1,Nx ),  ( beam%y(j), j=1,Ny )
      REWIND( unitNo )

      beam%dx = dx
      beam%dy = dy
      beam%wavelength = wavelength
      beam%curvature  = curvature

   ELSE
      WRITE( ioUnit, "(I3)" ) unitNo
      CALL ErrorExit( "RetrieveWavefront: error reading from unit "//ioUnit )
   END IF

   END SUBROUTINE RetrieveWavefront



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE WriteWavefrontUNF( beam, fileName, ioUnit )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Writes wavefront data to an unformatted file, for off-line analysis. The
   ! output file is given the name "fileName.unf" and it can be read by MATLAB.
   !
   !  USES NON-STANDARD FORTRAN 95
   !
   !		25 June '01
   !      7 Jan  '02    outputName may be 64 characters in length
   !     12 Mar  '02    Reversed order of writing array rows and columns
   !     10 June '02    Output file name up to 128 characters

   USE Utility_Routines

   IMPLICIT NONE

   TYPE(wavefront),   INTENT(IN) :: beam
   CHARACTER(LEN=*),  INTENT(IN) :: fileName
   INTEGER, OPTIONAL, INTENT(IN) :: ioUnit

   INTRINSIC SIZE, TRIM, PRESENT, REAL, AIMAG

   ! Local variables
   INTEGER            :: Nx, Ny
   CHARACTER(LEN=128) :: outputName
   INTEGER            :: unitNumber
   INTEGER, PARAMETER :: defaultUnit = 1
   INTEGER            :: ioFlag
   INTEGER, PARAMETER :: ok = 0
   INTEGER            :: j, k
   INTEGER            :: fileSize
   ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   Nx = SIZE( beam%amp, 2 )
   Ny = SIZE( beam%amp, 1 )

   ! Number of bytes for:
   !	- Nx and Ny ( 4 bytes each )
   !  - wavelength and curvature ( 8 bytes each )
   !  - x-coordinate vector ( 8*Nx )
   !  - y-coordinate vector ( 8*Ny )
   !  - complex beam 2-d array ( 2*8*Nx*Ny )
   !  => fileSize = 4*2 + 8*2 + 8*( Nx + Ny + 2*Nx*Ny )

   fileSize = 8*( 3 + Nx + Ny + 2*Nx*Ny )

   outputName = TRIM( fileName )//".unf"

   unitNumber = defaultUnit
   IF ( PRESENT( ioUnit ) ) unitNumber = ioUnit

   OPEN( unitNumber, FILE     = TRIM( outputName ),	&
						   STATUS   = "NEW",						&
						   FORM     = "BINARY",	  				&  ! NON-STANDARD
                     ACCESS   = "SEQUENTIAL",			&
                     BLOCKSIZE = fileSize,				&  ! NON-STANDARD
						   ACTION   = "WRITE",					&
						   POSITION = "REWIND",					&
						   IOSTAT   = ioFlag 					)

   IF ( ioFlag /= ok )  &
	   CALL ErrorExit( "WriteWavefrontUNF: error opening output file" )

   WRITE( unitNumber ) Nx, Ny
   WRITE( unitNumber ) beam%wavelength
   WRITE( unitNumber ) beam%curvature

   WRITE( unitNumber, IOSTAT=ioFlag ) (beam%x(k),k=1,Nx)
   IF ( ioFlag /= ok )  &
	   CALL ErrorExit( "WriteWavefrontUNF: error writing beam x values" )

   WRITE( unitNumber, IOSTAT=ioFlag ) (beam%y(j),j=1,Ny)
   IF ( ioFlag /= ok )  &
	   CALL ErrorExit( "WriteWavefrontUNF: error writing beam y values" )

   WRITE( unitNumber, IOSTAT=ioFlag ) ((REAL( beam%amp(j,k) ),j=1,Ny),k=1,Nx)
   IF ( ioFlag /= ok )  &
	   CALL ErrorExit( "WriteWavefrontUNF: error writing beam real part" )

   WRITE( unitNumber, IOSTAT=ioFlag ) ((AIMAG( beam%amp(j,k) ),j=1,Ny),k=1,Nx)
   IF ( ioFlag /= ok )  &
	   CALL ErrorExit( "WriteWavefrontUNF: error writing beam imaginary part" )

   CLOSE( unitNumber )

	END SUBROUTINE WriteWavefrontUNF

END MODULE Wavefronts
