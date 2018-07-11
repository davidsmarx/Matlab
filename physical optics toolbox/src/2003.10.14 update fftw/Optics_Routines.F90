!     Last change:  RSB  17 Mar 2003    2:07 pm
MODULE Optics_Routines

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   ! SIM-specific routines removed from this module 14 March 2003

	USE Kinds
   USE SI_Units
   USE Angle_Units
   USE Wavefronts
   USE Error_Exit

   IMPLICIT NONE

	! DEFINE LENS OBJECT....................................................
	TYPE, PUBLIC :: simpleLens
		PRIVATE
		REAL(pr) :: focLen
		REAL(pr) :: diam      = 0.0_pr
		REAL(pr) :: xDecenter = 0.0_pr
		REAL(pr) :: yDecenter = 0.0_pr
		TYPE (zernikeAberration) :: aberr
	END TYPE simpleLens
	!............................................................................

   PRIVATE SetLensSingleAberr, SetLensMultipleAberr

	INTERFACE SetLensAberration
		MODULE PROCEDURE SetLensSingleAberr, SetLensMultipleAberr
	END INTERFACE

CONTAINS


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Clip( beam, shape, size1, size2, xc, yc, angle, obsc, &
                    tiltAngle, tiltOrientation                       )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calls appropriate wavefront clipping routine, depending on shape of
   ! aperture or obscuration. Aperture may be tilted about an axis that makes
   ! a counter-clockwise angle with respect to the x-axis given by tiltOrientation
   !      5 Feb '02
   !		27 Feb '03		Added "annulus" clipping
   !     14 Mar '03     Aperture may be tilted

   !  beam     Wavefront to be clipped
   !  shape    Clipping shape--"circle" "square" "rectangle" or "annulus"
   !  size1    Diameter of circle, width of square or rectangle
   !  size2    OPTIONAL, height of aperture--needed for rectangle only
   !  xc, yc   OPTIONAL, location of aperture center
   !  angle    OPTIONAL, counterclockwise rotation of square or rectangle
   !  obsc     OPTIONAL, .TRUE. if obscuration--default is .FALSE.
   !  tiltAngle         OPTIONAL, amount of tilt
   !  tiltOrientation   OPTIONAL, orientation of tilt axis w.r.t. x-axis
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Kinds
   USE Wavefronts
   USE Error_Exit

   IMPLICIT NONE

   INTRINSIC PRESENT

   TYPE(wavefront),  INTENT(IN OUT)       :: beam
   CHARACTER(LEN=*), INTENT(IN)           :: shape
   REAL(pr),         INTENT(IN)           :: size1
   REAL(pr),         INTENT(IN), OPTIONAL :: size2
   REAL(pr),         INTENT(IN), OPTIONAL :: xc, yc
   REAL(pr),         INTENT(IN), OPTIONAL :: angle
   LOGICAL,          INTENT(IN), OPTIONAL :: obsc
   REAL(pr),         INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

   ! Local variables
   REAL(pr) :: xo, yo
   REAL(pr) :: rotAngle
   LOGICAL  :: obscuration
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	xo          = 0.0_pr
	yo          = 0.0_pr
	rotAngle    = 0.0_pr
   obscuration = .FALSE.
   IF ( PRESENT( xc ) )    xo          = xc
   IF ( PRESENT( yc ) )    yo          = yc
   IF ( PRESENT( angle ) ) rotAngle    = angle
   IF ( PRESENT( obsc ) )  obscuration = obsc

   SELECT CASE ( shape )
      CASE ( "circle", "Circle", "CIRCLE" )
         CALL ClipCirc( beam, size1, xo, yo, obscuration, &
                        tiltAngle, tiltOrientation         )
      CASE ( "square", "Square", "SQUARE" )
         CALL ClipRect( beam, size1, size1, xo, yo, rotAngle, obscuration, &
                        tiltAngle, tiltOrientation                          )
      CASE ( "rectangle", "Rectangle", "RECTANGLE" )
         IF ( .NOT.PRESENT( size2 ) )  &
            CALL ErrorExit( "Clip: missing rectangular aperture size parameter" )
         CALL ClipRect( beam, size1, size2, xo, yo, rotAngle, obscuration, &
                        tiltAngle, tiltOrientation                          )
      CASE ( "annulus", "Annulus", "ANNULUS" )
      	! Assumes concentric inner and outer limiting circles
         IF ( .NOT.PRESENT( size2 ) )  &
            CALL ErrorExit( "Clip: missing annular aperture size parameter" )
      	CALL ClipAnnulus( beam, size1, size2, xo, yo, xo, yo, obscuration, &
      	                  tiltAngle, tiltOrientation                        )
      CASE DEFAULT
         CALL ErrorExit( "Clip: unrecognized aperture shape" )
   END SELECT

   END SUBROUTINE Clip



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ClipAnnulus( beam, Douter, Dinner, xco, yco, xci, yci, obsc, &
                           tiltAngle, tiltOrientation                       )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips input wavefront by an annular aperture if obsc = .FALSE. (default);
   ! otherwise, clips by an annular obscuration. The outer and inner portions
   ! of the annulus are centered at (xco,yco) and (xci,yci), respectively;
   ! default centers are both (0,0). If (xco,yco) is specified but (xci,yci) is
   ! not, the inner portion is centered at (xco,yco). Tilt of the aperture is
   ! specified (optionally) by the amount of tilt (tiltAngle) and the orient-
   ! ation (tiltOrientation) of the tilt axis with respect to the x-axis.
   !  27 Feb '03
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Clipping_Routines
   USE Wavefronts

   IMPLICIT NONE

   INTRINSIC PRESENT, MAX

   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr),        INTENT(IN)     :: Douter, Dinner
   REAL(pr),        INTENT(IN), OPTIONAL :: xco, yco, xci, yci
   LOGICAL,         INTENT(IN), OPTIONAL :: obsc
   REAL(pr),        INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

   ! Local variables
   TYPE(wavefront) :: centerPart
   REAL(pr)        :: x1, x2, y1, y2
   REAL(pr)        :: dx, dy
   INTEGER         :: Nx, Ny
   REAL(pr)        :: xc, yc
   REAL(pr)        :: xcInner, ycInner
   REAL(pr)        :: halfWidth
   REAL(pr)        :: tilt, orientation
   REAL(pr)        :: radius
   LOGICAL         :: isObsc
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   isObsc = .FALSE.
   xc     = 0.0_pr
   yc     = 0.0_pr
   tilt        = 0.0_pr
   orientation = 0.0_pr

   IF ( PRESENT( obsc ) ) isObsc = obsc
   IF ( PRESENT( xco ) )      xc = xco
   IF ( PRESENT( yco ) )      yc = yco
   IF ( PRESENT( tiltAngle ) )       tilt        = tiltAngle
   IF ( PRESENT( tiltOrientation ) ) orientation = tiltOrientation
   xcInner = xc
   ycInner = yc
   IF ( PRESENT( xci ) ) xcInner = xci
   IF ( PRESENT( yci ) ) ycInner = yci

   IF ( isObsc ) THEN

      CALL GetWavefrontSampling( beam, Nx, Ny, dx, dy )
      CALL GetEdgeSmoothing( halfWidth )
      radius = 0.5_pr*Dinner + 2.0_pr*MAX( dx, dy )*( halfWidth + 1.0_pr )

      x1 = -radius + xcInner
      x2 =  radius + xcInner
      y1 = -radius + ycInner
      y2 =  radius + ycInner
      CALL GetWavefrontSection( beam, x1, x2, y1, y2, centerPart )

      CALL ClipCirc( beam, Douter, xc, yc, isObsc, tilt, orientation )
      CALL ClipCirc( centerPart, Dinner, xcInner, ycInner, .NOT.isObsc, &
                     tilt, orientation                                   )
      CALL AddWavefrontSection( beam, centerPart )
      CALL ClearWavefront( centerPart )

   ELSE

      CALL ClipCirc( beam, Douter, xc, yc, isObsc, tilt, orientation )
      CALL ClipCirc( beam, Dinner, xcInner, ycInner, .NOT.isObsc,  &
                     tilt, orientation                              )

   END IF

   END SUBROUTINE ClipAnnulus



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Lens( thisLens, beam, applyCurv )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Applies a lens (possibly decentered) to input wavefront. Clips by
   !  lens aperture (zero diameter value suppresses clipping) and adds Zernike
   !	aberrations, if any.
	!		 8 Oct '02
   !     11 Mar '03     Reversed order of beam and thisLens input arguments
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	IMPLICIT NONE

   TYPE(simpleLens),  INTENT(IN)   :: thisLens
   TYPE(wavefront), INTENT(IN OUT) :: beam
   LOGICAL, OPTIONAL, INTENT(IN)   :: applyCurv

   ! Local variables
   REAL(pr) :: xo, yo
   REAL(pr) :: focLen, diam
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   focLen = thisLens%focLen
   diam   = thisLens%diam
   xo     = thisLens%xDecenter
   yo     = thisLens%yDecenter

	CALL ThinLens( beam, focLen, diam, xo, yo, applyCurv )
   CALL ApplyAberration( beam, thisLens%aberr )

   END SUBROUTINE Lens


	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	 SUBROUTINE SetLensProperties( lens, focLen, diam, xDecenter, yDecenter )
	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	 ! Set focal length and diameter of a lens. Aberrations can be specified
	 ! by a call to SetLensAberr; default aberrations are zero.
	 !		 8 Oct '02
	 !

	IMPLICIT NONE

	TYPE( simpleLens ), INTENT(IN OUT) :: lens
	REAL(pr),           INTENT(IN)     :: focLen
	REAL(pr), OPTIONAL, INTENT(IN)     :: diam
	REAL(pr), OPTIONAL, INTENT(IN)     :: xDecenter, yDecenter

	INTRINSIC PRESENT, ABS, EPSILON

	! Local variables (none)
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( focLen ) > EPSILON( 0.0_pr ) ) THEN
		lens%focLen = focLen
   ELSE
      CALL ErrorExit( "SetThinLens: invalid focal length" )
   END IF

	IF ( PRESENT( diam ) )      lens%diam      = diam
	IF ( PRESENT( xDecenter ) ) lens%xDecenter = xDecenter
	IF ( PRESENT( yDecenter ) ) lens%yDecenter = yDecenter

	END SUBROUTINE SetLensProperties



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	SUBROUTINE SetLensMultipleAberr( lens, modes, coeffs, refDiam, RMS )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	! Set parameters describing thin lens aberrations.
	!		 7 Oct '02
	!
	! REQUIRED INPUT
	!
	! lens			A lens object
	!
	! modes()		Array of mode numbers for aberration
	! coeffs()		Array of Zernike coefficients (OPD values)
	!
	! OPTIONAL INPUT
	!
	! refDiam		Reference diameter for defining Zernike polynomials; if not
	!					present, physical diameter will be used
	! RMS				Input value .TRUE. if coefficients are RMS values
	!					(default is .FALSE., denoting peak-to-valley)

	IMPLICIT NONE

	TYPE( simpleLens ), INTENT(IN OUT) :: lens
	INTEGER,  DIMENSION(:), INTENT(IN) :: modes
	REAL(pr), DIMENSION(:), INTENT(IN) :: coeffs
	REAL(pr), OPTIONAL,     INTENT(IN) :: refDiam
	LOGICAL,  OPTIONAL,     INTENT(IN) :: RMS

	INTRINSIC PRESENT

	! Local variables:
	REAL(pr) :: xo, yo, rD
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	IF ( PRESENT( refDiam ) ) THEN
		rD = refDiam
	ELSE
		IF ( lens%diam > 0.0_pr ) THEN
			rD = lens%diam
		ELSE
			CALL ErrorExit( "SetThinLens: need Zernike reference diameter" )
		END IF
	END IF
	xo = lens%xDecenter
	yo = lens%yDecenter
	CALL SetZernikeAberration( lens%aberr, modes, coeffs, rD, xo, yo, RMS )

	END SUBROUTINE SetLensMultipleAberr



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	SUBROUTINE SetLensSingleAberr( lens, mode, coeff, refDiam, RMS )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	! Set parameters describing thin lens aberrations (a single Zernike mode)
	!		 7 Oct '02
	!
	! REQUIRED INPUT
	!
	! lens			A lens object
	!
	! mode			Mode number for this aberration
	! coeff			Zernike coefficient (OPD value)
	!
	! OPTIONAL INPUT
	!
	! refDiam		Reference diameter for defining Zernike polynomials; if not
	!					present, physical diameter will be used
	! RMS				.TRUE. if coefficient is RMS value
	!					(default is .FALSE., denoting peak-to-valley)

	IMPLICIT NONE

	TYPE( simpleLens ), INTENT(IN OUT) :: lens
	INTEGER,            INTENT(IN) :: mode
	REAL(pr),           INTENT(IN) :: coeff
	REAL(pr), OPTIONAL, INTENT(IN) :: refDiam
	LOGICAL,  OPTIONAL, INTENT(IN) :: RMS

	INTRINSIC PRESENT

	! Local variables:
	REAL(pr) :: xo, yo, rD
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	IF ( PRESENT( refDiam ) ) THEN
		rD = refDiam
	ELSE
		IF ( lens%diam > 0.0_pr ) THEN
			rD = lens%diam
		ELSE
			CALL ErrorExit( "SetThinLens: need Zernike reference diameter" )
		END IF
	END IF
	xo = lens%xDecenter
	yo = lens%yDecenter
	CALL SetZernikeAberration( lens%aberr, mode, coeff, rD, xo, yo, RMS )

	END SUBROUTINE SetLensSingleAberr



	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	 SUBROUTINE DecenterLens( lens, xDecenter, yDecenter )
	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	 ! Change decentering of input lens.
	 !		 9 Oct '02
	 !

	IMPLICIT NONE

	TYPE( simpleLens ), INTENT(IN OUT) :: lens
	REAL(pr),           INTENT(IN)     :: xDecenter, yDecenter

	! Local variables (none)
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	lens%xDecenter = xDecenter
	lens%yDecenter = yDecenter

	END SUBROUTINE DecenterLens



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE PropagateToFocalPlane( thisLens, beam,                &
   											 dxFp, dyFp, deFocus, applyCurv  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Applies curvature and WFE (if any) to wavefront, then propagates to the
   ! lens's focal plane.
   !      9 Oct  '02
   !     11 Mar  '03    Order of input beam and thisLens interchanged
   !
   ! REQUIRED INPUT
   !  thisLens   A lens object; should have positive focal length
   !  beam       Wavefront incident on lens; output is final-plane wavefront
   ! OPTIONAL INPUT
   !  dxFp, dyFp  Required sample spacing in wavefront grid in focal plane
   !  deFocus     Displacement of final plane from perfect focus
   !  applyCurv   Flag for applying curvature in focal plane (default .FALSE.)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(simpleLens),   INTENT(IN)  :: thisLens
   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr), OPTIONAL, INTENT(IN)  :: dxFp, dyFp
   REAL(pr), OPTIONAL, INTENT(IN)  :: deFocus
   LOGICAL,  OPTIONAL, INTENT(IN)  :: applyCurv

   INTRINSIC PRESENT

   ! Local variables:
   REAL(pr) :: dz
   LOGICAL  :: applCrv
   REAL(pr) :: focLen
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   focLen = thisLens%focLen
   IF ( focLen < 0 )  &
   	CALL ErrorExit( "PropagateToFocalPlane: negative focal length?" )
   CALL Lens( thisLens, beam )

   dz      = 0.0_pr
   applCrv = .FALSE.
   IF ( PRESENT( deFocus ) )   dz      = deFocus
   IF ( PRESENT( applyCurv ) ) applCrv = applyCurv

   IF ( PRESENT( dxFp )  .AND.  PRESENT( dyFp ) ) THEN
      CALL Propagate( beam, focLen+dz, dxFp, dyFp, applCrv )
   ELSE
      ! Use default (i.e., FFT) sample spacing in far field
      CALL Propagate( beam, focLen+dz, applyCurv=applCrv )
   END IF

   END SUBROUTINE PropagateToFocalPlane



	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    FUNCTION LensFocalLength( lens ) RESULT( focLen )
	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    ! Returns value of focal length for a lens object
    !     9 Oct  '02
    !
	 !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(simpleLens), INTENT(IN) :: lens
   REAL(pr)                     :: focLen

	! Local variables (none)
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   focLen = lens%focLen

   END FUNCTION LensFocalLength



	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    FUNCTION LensDiameter( lens ) RESULT( diam )
	!= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

    ! Returns value of diameter for a lens object
    !     9 Oct  '02
    !
	 !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(simpleLens), INTENT(IN) :: lens
   REAL(pr)                     :: diam

	! Local variables (none)
	!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   diam = lens%diam

   END FUNCTION LensDiameter



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ThinLens( beam, focalLength,                  &
   							diameter, xDecen, yDecen, applyCurv  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Applies thin lens (possibly decentered) to input wavefront. Clips by
   !  lens aperture (zero diameter value suppresses clipping). No aberrations.
   !  Called by Lens subroutine (9 Oct '02).
	!		11 Jun  '01
   !      5 July '01    Will apply residual curvature according to "applyCurv",
   !							but does not apply curvature if applyCurv not present
   !     19 Dec  '01    Diameter is optional input; other minor changes
   !		 7 Oct  '02		Corrected logic for applying piston and tilt.
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	IMPLICIT NONE

   TYPE(wavefront), INTENT(IN OUT) :: beam
   REAL(pr),            INTENT(IN) :: focalLength
   REAL(pr), OPTIONAL,  INTENT(IN) :: diameter
   REAL(pr), OPTIONAL,  INTENT(IN) :: xDecen, yDecen
   LOGICAL,  OPTIONAL,  INTENT(IN) :: applyCurv

   INTRINSIC PRESENT

   ! Local variables
   REAL(pr) :: xo,   yo
   REAL(pr) :: curv, xTilt, yTilt, piston
   LOGICAL  :: fullCurvature
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Positive focal length produces negative wavefront curvature.
   ! The beam's wavefront curvature variable is modified, but the curvature
   ! is not explicitly applied to the complex field values, unless "applyCurv"
   ! is input as .TRUE.

   curv = -1.0_pr/focalLength
   CALL ModifyCurvature( beam, curv )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xDecen ) ) xo = xDecen
   IF ( PRESENT( yDecen ) ) yo = yDecen

   IF ( PRESENT( applyCurv ) ) THEN
   	fullCurvature = applyCurv
   ELSE
   	fullCurvature = .FALSE.
   END IF

   IF ( fullCurvature ) THEN
		CALL ApplyCurvature( beam, xc=xo, yc=yo )
   ELSE IF ( ABS( xo ) > EPSILON( 0.0_pr )  .OR.  &
   			 ABS( yo ) > EPSILON( 0.0_pr )        ) THEN
   	! A decentered lens introduces piston and tilt, which are applied
   	! to the complex field when the full curvature is not applied.
   	xTilt  = -xo*curv
   	yTilt  = -yo*curv
   	piston = -0.5_pr*(xo*xTilt + yo*yTilt)
		curv   =  0.0_pr
		CALL PistonTiltFocus( beam, piston, xTilt, yTilt, curv )
   END IF

   IF ( PRESENT( diameter ) ) CALL ClipCirc( beam, diameter, xo, yo )

   END SUBROUTINE ThinLens



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE HoleGrating( beam, pattern, holeSpacing, holeShape, size1,  &
   								size2, xc, yc, angle                            )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips input wavefront using a regular array of apertures.
   !      8 Apr '02
   !     22 Apr '02     Grating may be rotated (optional input)

   !  beam     	Wavefront to be clipped
   !  pattern  	Arrangement of apertures: "square" or "hexagonal"
   !  holeSpacing Distance between aperture centers (nearest neighbor)
   !  holeShape   Aperture shape--"circle", "square" or "rectangle"
   !  size1    	Diameter of circle, width of square or rectangle
   !  size2    	OPTIONAL, height of aperture--needed for rectangle only
   !  xc, yc   	OPTIONAL, location of array center (default 0.0,0.0)
   !  angle       OPTIONAL, counter-clockwise rotation
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Kinds
   USE Constants
   USE Array_Routines,    ONLY : CoordinateVector, R_2d
   USE Clipping_Routines, ONLY : GetEdgeSmoothing
   USE Polygons
   USE Wavefronts
   USE Error_Exit

   IMPLICIT NONE

   INTRINSIC FLOOR, PRESENT

   TYPE(wavefront), INTENT(IN OUT) :: beam

   CHARACTER(LEN=*) :: pattern
   REAL(pr)         :: holeSpacing
   CHARACTER(LEN=*) :: holeShape
   REAL(pr)         :: size1

   REAL(pr), OPTIONAL :: size2
   REAL(pr), OPTIONAL :: xc, yc
   REAL(pr), OPTIONAL :: angle

   ! Local variables
   LOGICAL  :: isHex
   INTEGER  :: NxHoles, kx
   INTEGER  :: NyHoles, jy
   REAL(pr) :: xSize,   ySize
   REAL(pr) :: xSpacing, ySpacing
   REAL(pr) :: xo,  yo
   REAL(pr) :: theta
   INTEGER  :: Nxb, Nyb
   REAL(pr) :: dxb, dyb
   REAL(pr) :: xbo, ybo
   REAL(pr) :: halfWidth
   REAL(pr) :: xWmin, xWmax, yWmin, yWmax
   REAL(pr) :: xMin,  xMax,  yMin,  yMax
   TYPE(polygon) :: grid
   TYPE(polygon) :: rect
   REAL(pr), DIMENSION(2,2) :: gridBounds
   REAL(pr), DIMENSION(2,2) :: rectBounds
   REAL(pr), DIMENSION(2,2) :: rotMatrix
   REAL(pr), DIMENSION(2)   :: holeCenter
   REAL(pr), DIMENSION(:), ALLOCATABLE :: xHole, yHole
   REAL(pr), DIMENSION(:), ALLOCATABLE :: xGrid, yGrid
   TYPE(wavefront)         :: tempBeam
   TYPE(wavefront)         :: section
   COMPLEX(prc), PARAMETER :: cmplxZero = ( 0.0_pr, 0.0_pr )
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xSpacing = holeSpacing

   SELECT CASE ( pattern )
      CASE ( "hexagonal", "Hexagonal", "HEXAGONAL" )
         isHex = .TRUE.
         ySpacing = 0.5_pr*SQRT3*holeSpacing
      CASE ( "square",    "Square",    "SQUARE" )
         isHex = .FALSE.
		   ySpacing = holeSpacing
      CASE DEFAULT
         CALL ErrorExit( "HoleGrating: incorrect array shape parameter" )
   END SELECT

   xSize = size1

   SELECT CASE ( holeShape )
      CASE ( "circle",    "Circle",    "CIRCLE" )
   		ySize = size1
      CASE ( "square",    "Square",    "SQUARE" )
		   ySize = size1
      CASE ( "rectangle", "Rectangle", "RECTANGLE" )
         IF ( .NOT.PRESENT( size2 ) )  &
            CALL ErrorExit( "HoleGrating: rectangle dimension not specified" )
         ySize = size2
      CASE DEFAULT
         CALL ErrorExit( "HoleGrating: incorrect hole shape parameter" )
   END SELECT

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   theta = 0.0_pr
   IF ( PRESENT( angle ) ) theta = angle

   CALL GetWavefrontSampling( beam, Nxb, Nyb, dxb, dyb, xbo, ybo )
   ALLOCATE( xGrid(Nxb), yGrid(Nyb) )
   xGrid = CoordinateVector( Nxb, dxb, xbo )
   yGrid = CoordinateVector( Nyb, dyb, ybo )

   ! Determine height and width of beam grid, rotated by -theta
   grid  = MakeRectangle( xGrid(Nxb)-xGrid(1), yGrid(Nyb)-yGrid(1), -theta )
   gridBounds = BoundingRectangle( grid )
   CALL ClearPolygon( grid )

   ! Array of holes will be large enough to fill the final grid
   NxHoles = 2*FLOOR( (gridBounds(1,2)-gridBounds(1,1))/(2.0_pr*xSpacing) ) + 1
   NyHoles = 2*FLOOR( (gridBounds(2,2)-gridBounds(2,1))/(2.0_pr*ySpacing) ) + 1
   ALLOCATE( xHole(NxHoles), yHole(NyHoles) )

   ! These are positions of the holes before the array is rotated
   xHole = CoordinateVector( NxHoles, xSpacing, xo )
   yHole = CoordinateVector( NyHoles, ySpacing, yo )


   ! Find the rectangle large enough to enclose a single hole, including a
   ! small guard band. This will be used to get a section of the initial
   ! wavefront for clipping by an individual hole.
   rect = MakeRectangle( xSize, ySize, theta )
   rectBounds = BoundingRectangle( rect )
   CALL ClearPolygon( rect )

   CALL GetEdgeSmoothing( halfWidth )
   halfWidth = halfWidth + 1.1_pr            ! Add a little extra...
   xWmin = rectBounds(1,1) - halfWidth*dxb
   xWmax = rectBounds(1,2) + halfWidth*dxb
   yWmin = rectBounds(2,1) - halfWidth*dyb
   yWmax = rectBounds(2,2) + halfWidth*dyb

   ! For each of the holes, select the piece of original wavefront that is
   ! transmitted,and add it to the new wavefront, whose amplitude is initially
   ! zero

   tempBeam = beam
   CALL SetWavefrontAmplitude( beam, cmplxZero )
   rotMatrix = R_2d( theta )

   DO kx = 1,NxHoles
      DO jy = 1,NyHoles,2     			! Do odd-numbered rows

         ! Calculate this hole's position after the array is rotated
         holeCenter = MATMUL( RotMatrix, (/ xHole(kx), yHole(jy) /) )

	      xMin = holeCenter(1) + xWmin
	      xMax = holeCenter(1) + xWmax
         yMin = holeCenter(2) + yWmin
         yMax = holeCenter(2) + yWmax

         ! Skip to next hole if this one is too close to edge of beam grid
         IF ( xMin < xGrid(1)  .OR.  xMax > xGrid(Nxb)  .OR. &
              yMin < yGrid(1)  .OR.  yMax > yGrid(Nyb)       ) CYCLE

         CALL GetWavefrontSection( tempBeam, xMin, xMax, yMin, yMax, section )
         CALL Clip( section, holeShape, xSize, ySize,   &
         			  holeCenter(1), holeCenter(2), theta )
         CALL AddWavefrontSection( beam, section )

      END DO
   END DO

   IF ( isHex ) THEN                   ! For a hexagonal array, there is
      NxHoles = NxHoles - 1            ! one less hole in the even rows,
      xHole = xHole + 0.5_pr*xSpacing  ! and centers are shifted
   END IF

   DO kx = 1,NxHoles
      DO jy = 2,NyHoles-1,2				! Do even-numbered rows

         holeCenter = MATMUL( RotMatrix, (/ xHole(kx), yHole(jy) /) )

	      xMin = holeCenter(1) + xWmin
	      xMax = holeCenter(1) + xWmax
         yMin = holeCenter(2) + yWmin
         yMax = holeCenter(2) + yWmax

         IF ( xMin < xGrid(1)  .OR.  xMax > xGrid(Nxb)  .OR. &
              yMin < yGrid(1)  .OR.  yMax > yGrid(Nyb)       ) CYCLE

         CALL GetWavefrontSection( tempBeam, xMin, xMax, yMin, yMax, section )
         CALL Clip( section, holeShape, xSize, ySize,   &
         			  holeCenter(1), holeCenter(2), theta )
         CALL AddWavefrontSection( beam, section )

      END DO
   END DO

   CALL ClearWavefront( tempBeam )
   CALL ClearWavefront( section )
   DEALLOCATE( xHole, yHole )
   DEALLOCATE( xGrid, yGrid )

   END SUBROUTINE HoleGrating

END MODULE Optics_Routines
