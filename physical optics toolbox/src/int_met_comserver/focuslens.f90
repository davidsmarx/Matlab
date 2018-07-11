MODULE ModFocusLens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Propagates to lens, Applies lens and propagates one focal length behind lens


   USE Kinds
   USE SI_Units
   USE Constants
   USE Wavefronts
   USE Optics_Routines
   USE Utility_Routines


   IMPLICIT NONE

   PRIVATE


    !! AVAILABLE SUBROUTINES:
    PUBLIC :: FocusLens


    !!!!!!!!!!!!!!!!!!!!!!!! Default Internal Metrology System Parameters !!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!! Basic Beam Parameters
	INTEGER, PARAMETER  :: Nx = 1024
	INTEGER, PARAMETER  :: Ny = 1024                  ! grid size
	REAL(pr), PARAMETER :: Gx = 25_pr*mm
	REAL(pr), PARAMETER :: Gy = 25_pr*mm                ! beam extent
	REAL(pr), PARAMETER :: lambda = 1.3_pr*micrometers  ! Beam wavelength
	

    !! Heterodyne section parameters.
	REAL(pr), PARAMETER :: fieldsep_z = 35_pr*mm     ! meters distance to next component (focusing lens)
	
	!! Focusing lens
	REAL(pr), PARAMETER :: hetfocuslens_f = 35_pr*mm    ! focal length
	REAL(pr), PARAMETER :: hetfocuslens_D = 15_pr*mm    ! clear aperture of lens
	REAL(pr), PARAMETER :: hetfocuslens_z = hetfocuslens_f  ! distance to detector

	!! Field Integrating Detector
	CHARACTER(LEN=6)    :: hetdetector_shape = "circle"
	REAL(pr), PARAMETER :: hetdetector_D = 0.3_pr*mm    ! detector diameter

   !! unit #'s for outputs
   INTEGER, PARAMETER :: iUNF        = 1   ! ionumber for writing UNF outputs

    ! Sample spacing is derived from grid size
!	REAL(pr) :: dx = Gx/REAL( Nx, KIND=pr ) 
!	REAL(pr) :: dy = Gy/REAL( Ny, KIND=pr )

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FocusLens( beam )

IMPLICIT NONE

   type(Wavefront), INTENT(INOUT) :: beam  ! the beam

   TYPE(simpleLens)    :: hetfocuslens              ! lens object
   REAL(pr)            :: dxout, dyout, dx, dy
   INTEGER             :: Nx, Ny

   CALL GetWavefrontSampling( beam, Nx, Ny, dx, dy )


    !! propagate to the next element (focusing lens)
    call Propagate( beam, fieldsep_z, dx, dy, applyCurv=.TRUE. )

   !! initialize focus lens and apply focusing lens as a subroutine
   call SetLensProperties( hetfocuslens, hetfocuslens_f, hetfocuslens_D )
   dxout = 1.05_pr*hetdetector_D/REAL( Nx, KIND = pr )
   dyout = 1.05_pr*hetdetector_D/REAL( Ny, KIND = pr )

   call PropagateToFocalPlane( hetfocuslens, beam, dxout, dyout, applyCurv=.TRUE. )


END SUBROUTINE FocusLens

END MODULE ModFocusLens
