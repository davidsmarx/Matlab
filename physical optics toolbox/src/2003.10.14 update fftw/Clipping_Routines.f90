!     Last change:  RSB  14 Oct 2003   11:29 am
MODULE Clipping_Routines

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds
   USE Error_Exit

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GaussianCircularApodization
   PUBLIC :: CosineCircularApodization
   PUBLIC :: QuadraticCircularApodization
   PUBLIC :: SetEdgeSmoothing, GetEdgeSmoothing,  &
             RestoreDefaultEdgeSmoothing
   PUBLIC :: CircularClip,   CircularObsc,   &
   			 EllipticalClip, EllipticalObsc, &
             PolygonClip,    PolygonObsc

!!!!!!!! Additional Routines not in the official release
   PUBLIC :: WindowClip
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   INTERFACE CircularClip
      MODULE PROCEDURE CircularClipReal,   &
                       CircularClipCmplx
   END INTERFACE

   INTERFACE CircularObsc
      MODULE PROCEDURE CircularObscReal,   &
                       CircularObscCmplx
   END INTERFACE

   INTERFACE EllipticalClip
      MODULE PROCEDURE EllipticalClipReal,   &
                       EllipticalClipCmplx
   END INTERFACE

   INTERFACE EllipticalObsc
      MODULE PROCEDURE EllipticalObscReal,   &
                       EllipticalObscCmplx
   END INTERFACE

   INTERFACE PolygonClip
      MODULE PROCEDURE PolygonClipReal,   &
                       PolygonClipCmplx
   END INTERFACE

   INTERFACE PolygonObsc
      MODULE PROCEDURE PolygonObscReal,   &
                       PolygonObscCmplx
   END INTERFACE

   INTERFACE SoftClip
      MODULE PROCEDURE SoftClipReal,   &
                       SoftClipCmplx
   END INTERFACE

   ! Parameter for soft-clipping at aperture edges
   REAL(pr), PARAMETER :: defaultSmoothing = 1.0_pr   ! Default value
   REAL(pr) :: SMOOTHING_HALF_WIDTH = defaultSmoothing

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GaussianCircularApodization( array, x, y,                 &
                                           outerDiam, innerDiam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Apodizes array by Gaussian taper outside innerDiameter centered at (xc,yc).
   ! The amplitude falls to 1/e at radius equal to outerDiam/2.
   !     11 Mar '02
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x,  y
   REAL(pr),                     INTENT(IN)     :: outerDiam
   REAL(pr),                     INTENT(IN), OPTIONAL :: innerDiam
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC PRESENT, SQRT, SIZE, REAL, AIMAG, CMPLX

   ! Local variables
   REAL(pr) :: xo, yo
   REAL(pr) :: width, r_inner, xsq, r
   INTEGER  :: j, k
   REAL(pr) :: factor
   REAL(pr) :: reA, imA
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   IF ( PRESENT( innerDiam ) ) THEN
      r_inner = 0.5_pr*innerDiam
   ELSE
      r_inner = 0.0_pr
   END IF
   width = 0.5_pr*outerDiam - r_inner

   DO k = 1,SIZE( x )
      xsq = ( x(k) - xo )**2
   	DO j = 1,SIZE( y )
         r = SQRT( ( y(j)-yo )**2 + xsq )
         IF ( r > r_inner ) THEN
            factor = exp( -(( r - r_inner )/width)**2 )
            reA    = factor * REAL( array(j,k) )
            imA    = factor * AIMAG( array(j,k) )
            array(j,k) = CMPLX( reA, imA, KIND = prc )
         END IF
      END DO
   END DO

   END SUBROUTINE GaussianCircularApodization



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE CosineCircularApodization( array, x, y,                 &
                                         outerDiam, innerDiam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Apodizes array by cosine taper, varying from 1 to 0 outside innerDiameter
   ! centered at (xc,yc). The amplitude falls to 1/2 at radius equal to
   ! (outerDiam + innerDiam)/4, and is zero outside outerDiam.
   !     11 Mar '02
   !		18 Apr '02		Corrected: amplitude goes to zero at outerDiam; no
   !                    apodization if outerDiam is not greater than innerDiam.
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x,  y
   REAL(pr),                     INTENT(IN)     :: outerDiam
   REAL(pr),                     INTENT(IN)     :: innerDiam
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC PRESENT, EPSILON, SQRT, SIZE, REAL, AIMAG, CMPLX

   ! Local variables
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )
   REAL(pr) :: xo, yo
   REAL(pr) :: width, r_inner, r_outer, xsq, r
   INTEGER  :: j, k
   REAL(pr) :: factor
   REAL(pr) :: reA, imA
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( outerDiam < 0.0_pr  .OR.  innerDiam < 0.0_pr ) THEN
      CALL ErrorExit( "Negative diameter value(s) in Circ_Cos_Apodization" )
   ELSE IF ( (outerDiam - innerDiam) <= innerDiam*EPSILON( outerDiam ) ) THEN
      RETURN
   END IF

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   r_inner = 0.5_pr*innerDiam
   r_outer = 0.5_pr*outerDiam
   width   = r_outer - r_inner

   DO k = 1,SIZE( x )
      xsq = ( x(k) - xo )**2
   	DO j = 1,SIZE( y )
         r = SQRT( ( y(j)-yo )**2 + xsq )
         IF ( r >= r_outer ) THEN
            array(j,k) = ZERO
         ELSE IF ( r > r_inner ) THEN
            factor = cos( HALFPI*( r - r_inner )/width )**2
            reA    = factor * REAL( array(j,k) )
            imA    = factor * AIMAG( array(j,k) )
            array(j,k) = CMPLX( reA, imA, KIND = prc )
         END IF
      END DO
   END DO

   END SUBROUTINE CosineCircularApodization



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE QuadraticCircularApodization( array, x, y,                 &
   													  outerDiam, innerDiam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Apodizes array by the quadratic taper that is used for soft-clipping
   ! apertures. The apodization region is centered at (xc,yc). The amplitude
   ! falls to 1/2 at radius equal to (outerDiam + innerDiam)/4, and is zero
   ! outside outerDiam.
   !     10 Oct '02
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x,  y
   REAL(pr),                     INTENT(IN)     :: outerDiam
   REAL(pr),                     INTENT(IN)     :: innerDiam
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC PRESENT, EPSILON, SQRT, SIZE

   ! Local variables
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )
   REAL(pr) :: xo, yo
   REAL(pr) :: width, r_inner, r_outer, r_mid, xsq, r
   INTEGER  :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( outerDiam < 0.0_pr  .OR.  innerDiam < 0.0_pr ) THEN
      CALL ErrorExit( "Negative diameter value(s) in Quad_Circ_Apodization" )
   ELSE IF ( (outerDiam - innerDiam) <= innerDiam*EPSILON( outerDiam ) ) THEN
      RETURN
   END IF

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   r_inner = 0.5_pr*innerDiam
   r_outer = 0.5_pr*outerDiam
   r_mid   = 0.5_pr*( r_inner + r_outer )
   width   = r_outer - r_inner

   DO k = 1,SIZE( x )
      xsq = ( x(k) - xo )**2
   	DO j = 1,SIZE( y )
         r = SQRT( ( y(j)-yo )**2 + xsq )
         IF ( r >= r_outer ) THEN
            array(j,k) = ZERO
         ELSE IF ( r > r_inner ) THEN
            array(j,k) = SoftClip( array(j,k), ( r-r_mid )/width )
         END IF
      END DO
   END DO

   END SUBROUTINE QuadraticCircularApodization



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SetEdgeSmoothing( halfWidth )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Assigns a value to the smoothing half-width used by clipping routines.
   ! Its value probably should lie between 1.0 and 2.0.
   !     19 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr) :: halfWidth
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( halfWidth < 0.0_pr ) &
      CALL ErrorExit( "SetEdgeSmoothing: unphysical input value" )

   SMOOTHING_HALF_WIDTH = halfWidth

   END SUBROUTINE SetEdgeSmoothing



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE GetEdgeSmoothing( halfWidth )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns the current value of smoothing half-width.
   !     11 Feb  '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr) :: halfWidth
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   halfWidth = SMOOTHING_HALF_WIDTH

   END SUBROUTINE GetEdgeSmoothing



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE RestoreDefaultEdgeSmoothing
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Sets edge smoothing parameter to its default value
   !      8 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SMOOTHING_HALF_WIDTH = defaultSmoothing

   END SUBROUTINE RestoreDefaultEdgeSmoothing



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE CircularClipReal( array, x, y, diam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips real-valued array by circular aperture centered at (xc,yc).
   !     26 Oct  '00
   !      8 Feb  '02    Tests for aperture being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr), DIMENSION(:),   INTENT(IN)     :: x, y
   REAL(pr),                 INTENT(IN)     :: diam
   REAL(pr),                 INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, PRESENT

   ! Local variables
   REAL(pr), PARAMETER :: ZERO = 0.0_pr
   REAL(pr) :: dx, dy
   REAL(pr) :: radius, dr, r_inner, r_outer, xsq, r
   REAL(pr) :: xo, yo
   INTEGER  :: Ny, Nx, j, k
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, diam, xo, yo ) ) RETURN

   radius  = 0.5_pr * diam
   dr      = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   r_inner = radius - dr
   r_outer = radius + dr

   DO k = 1,Nx
      xsq = ( x(k) - xo )**2
   	DO j = 1,Ny
         r = SQRT( ( y(j)-yo )**2 + xsq )
         IF ( r >= r_outer ) THEN
            array(j,k) = ZERO
         ELSE IF ( r > r_inner ) THEN
            array(j,k) = SoftClip( array(j,k), ( r - radius )/dr )
         END IF
      END DO
   END DO

   END SUBROUTINE CircularClipReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE CircularClipCmplx( array, x, y, diam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips array by circular aperture centered at (xc,yc).
   !     26 Oct  '00
   !      8 Feb  '02    Tests for aperture being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !     26 Sept '02    Error corrected in argument of SoftClip
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x,  y
   REAL(pr),                     INTENT(IN)     :: diam
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, REAL, AIMAG, CMPLX

   ! Local variables
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )
   REAL(pr) :: dx, dy
   REAL(pr) :: xo, yo
   REAL(pr) :: radius, dr, r_inner, r_outer, xsq, r
   INTEGER  :: Ny, Nx, j, k
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, diam, xo, yo ) ) RETURN

   radius  = 0.5_pr * diam
   dr      = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   r_inner = radius - dr
   r_outer = radius + dr

   DO k = 1,Nx
      xsq = ( x(k) - xo )**2
   	DO j = 1,Ny
         r = SQRT( ( y(j)-yo )**2 + xsq )
         IF ( r >= r_outer ) THEN
            array(j,k) = ZERO
         ELSE IF ( r > r_inner ) THEN
            array(j,k) = SoftClip( array(j,k), ( r - radius )/dr )
         END IF
      END DO
   END DO

   END SUBROUTINE CircularClipCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE CircularObscReal( array, x, y, diam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips real-valued array by circular obscuration centered at (xc,yc).
   !     26 Oct  '00
   !      8 Feb  '02    Tests for obscuration being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr), DIMENSION(:),   INTENT(IN)     :: x,  y
   REAL(pr),                 INTENT(IN)     :: diam
   REAL(pr),                 INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, PRESENT

   ! Local variables
   REAL(pr), PARAMETER :: ZERO = 0.0_pr
   REAL(pr) :: dx, dy
   REAL(pr) :: xo, yo
   REAL(pr) :: radius, dr, r_inner, r_outer, xsq, r
   INTEGER  :: Nx, Ny, j, k
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, diam, xo, yo ) )  &
      CALL ErrorExit( "CircularObscReal: beam is completely obscured" )

   radius  = 0.5_pr * diam
   dr      = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   r_inner = radius - dr
   r_outer = radius + dr

   DO k = 1,Nx
      xsq = ( x(k) - xo )**2
   	DO j = 1,Ny
         r = SQRT( ( y(j)-yo )**2 + xsq )
         IF ( r <= r_inner ) THEN
            array(j,k) = ZERO
         ELSE IF ( r < r_outer ) THEN
            array(j,k) = SoftClip( array(j,k), ( radius - r )/dr )
         END IF
      END DO
   END DO

   END SUBROUTINE CircularObscReal


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE CircularObscCmplx( array, x, y, diam, xc, yc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips array by circular obscuration centered at (xc,yc).
   !     26 Oct  '00
   !      8 Feb  '02    Tests for obscuration being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x,  y
   REAL(pr),                     INTENT(IN)     :: diam
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, PRESENT, REAL, AIMAG, CMPLX

   ! Local variables
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )
   REAL(pr) :: dx, dy
   REAL(pr) :: xo, yo
   REAL(pr) :: radius, dr, r_inner, r_outer, xsq, r
   INTEGER  :: Nx, Ny, j, k
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, diam, xo, yo ) )  &
      CALL ErrorExit( "CircularObscComplex: beam is completely obscured" )

   radius  = 0.5_pr * diam
   dr      = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   r_inner = radius - dr
   r_outer = radius + dr

   DO k = 1,Nx
      xsq = ( x(k) - xo )**2
   	DO j = 1,Ny
         r = SQRT( ( y(j)-yo )**2 + xsq )
         IF ( r <= r_inner ) THEN
            array(j,k) = ZERO
         ELSE IF ( r < r_outer ) THEN
            array(j,k) = SoftClip( array(j,k), ( radius - r )/dr )
         END IF
      END DO
   END DO

   END SUBROUTINE CircularObscCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE EllipticalClipReal( array, x, y, majorDiam, minorDiam,  &
											 orientation, xc, yc                 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips real-valued array by elliptical aperture centered at (xc,yc).
   ! majorDiam		Major axis of ellipse
   ! minorDiam		Minor axis
   ! orientation	Angle (counter-clockwise) between grid x-axis
   !					and major axis of ellipse
   !     16 Sept '02		Based on circular clip routine
   !     22 Nov  '02    Calls revised EllipticalIteration
   !     23 Jul '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : R_2d

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr), DIMENSION(:),   INTENT(IN)     :: x, y
   REAL(pr),                 INTENT(IN)     :: majorDiam, minorDiam
   REAL(pr),                 INTENT(IN)     :: orientation
   REAL(pr),                 INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, PRESENT, MATMUL, SIGN

   ! Local variables
   REAL(pr), PARAMETER :: ZERO = 0.0_pr
   REAL(pr) :: dx, dy
   REAL(pr) :: dr
   REAL(pr) :: a, a_inner, a_outer
   REAL(pr) :: b, b_inner, b_outer
   REAL(pr) :: baSqr
   REAL(pr) :: xo, yo
   REAL(pr) :: xe, ye
   REAL(pr) :: edgeDist
   INTEGER  :: Ny, Nx, j, k
   REAL(pr), DIMENSION(2)   :: xy
   REAL(pr), DIMENSION(2,2) :: rotMatrix
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, minorDiam, xo, yo ) ) RETURN

   dr = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   a  = 0.5_pr * majorDiam
   b  = 0.5_pr * minorDiam
   a_inner = a - 2.0_pr*dr
   a_outer = a + 2.0_pr*dr
   b_inner = b - 2.0_pr*dr
   b_outer = b + 2.0_pr*dr

   baSqr = (b/a)**2

	rotMatrix = R_2d( -orientation )	! Transforms (x,y) -> (x',y'), a frame
												! in which x' and y' are along the major
												! and minor axes of the ellipse

   DO k = 1,Nx
   	DO j = 1,Ny
   		xy = MATMUL( rotMatrix, (/ x(k)-xo, y(j)-yo /) )

         IF ( (xy(1)/a_outer)**2 + (xy(2)/b_outer)**2 > 1.0_pr ) THEN
				! This point is completely outside the aperture
            array(j,k) = ZERO
         ELSE IF ( (xy(1)/a_inner)**2 + (xy(2)/b_inner)**2 > 1.0_pr ) THEN
         	! This point is in the "fringe" region; do edge smoothing
         	CALL EllipticalIteration( xy(1), xy(2), a, b, xe, ye )
				edgeDist = SQRT( (xy(1)-xe)**2 + (xy(2)-ye)**2 )
				edgeDist = SIGN( edgeDist, (xy(1)/a)**2 + (xy(2)/b)**2 - 1.0_pr )
				array(j,k) = SoftClip( array(j,k), edgeDist/dr )
         END IF

      END DO
   END DO

   END SUBROUTINE EllipticalClipReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE EllipticalClipCmplx( array, x, y, majorDiam, minorDiam,  &
											  orientation, xc, yc                 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips complex-valued array by elliptical aperture centered at (xc,yc).
   ! majorDiam		Major axis of ellipse
   ! minorDiam		Minor axis
   ! orientation	Angle (counter-clockwise) between grid x-axis
   !					and major axis of ellipse
   !     17 Sept '02		Based on circular clip routine
   !     22 Nov  '02    Calls revised EllipticalIteration
   !     23 July '03    Either xc or yc (or both) may be omitted from input
   !     14 Oct  '03    ZERO declared COMPLEX (rather than REAL)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : R_2d

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x, y
   REAL(pr),                     INTENT(IN)     :: majorDiam, minorDiam
   REAL(pr),                     INTENT(IN)     :: orientation
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, PRESENT, MATMUL, SIGN

   ! Local variables
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )
   REAL(pr) :: dx, dy
   REAL(pr) :: dr
   REAL(pr) :: a, a_inner, a_outer
   REAL(pr) :: b, b_inner, b_outer
   REAL(pr) :: baSqr
   REAL(pr) :: xo, yo
   REAL(pr) :: xe, ye
   REAL(pr) :: edgeDist
   INTEGER  :: Ny, Nx, j, k
   REAL(pr), DIMENSION(2)   :: xy
   REAL(pr), DIMENSION(2,2) :: rotMatrix
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, minorDiam, xo, yo ) ) RETURN

   dr = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   a  = 0.5_pr * majorDiam
   b  = 0.5_pr * minorDiam
   a_inner = a - 2.0_pr*dr
   a_outer = a + 2.0_pr*dr
   b_inner = b - 2.0_pr*dr
   b_outer = b + 2.0_pr*dr

   baSqr = (b/a)**2

	rotMatrix = R_2d( -orientation )	! Transforms (x,y) -> (x',y'), a frame
												! in which x' and y' are along the major
												! and minor axes of the ellipse

   DO k = 1,Nx
   	DO j = 1,Ny
   		xy = MATMUL( rotMatrix, (/ x(k)-xo, y(j)-yo /) )

         IF ( (xy(1)/a_outer)**2 + (xy(2)/b_outer)**2 > 1.0_pr ) THEN
				! This point is completely outside the aperture
            array(j,k) = ZERO
         ELSE IF ( (xy(1)/a_inner)**2 + (xy(2)/b_inner)**2 > 1.0_pr ) THEN
         	! This point is in the "fringe" region; do edge smoothing
         	CALL EllipticalIteration( xy(1), xy(2), a, b, xe, ye )
				edgeDist = SQRT( (xy(1)-xe)**2 + (xy(2)-ye)**2 )
				edgeDist = SIGN( edgeDist, (xy(1)/a)**2 + (xy(2)/b)**2 - 1.0_pr )
				array(j,k) = SoftClip( array(j,k), edgeDist/dr )
         END IF

      END DO
   END DO

   END SUBROUTINE EllipticalClipCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE EllipticalObscReal( array, x, y, majorDiam, minorDiam,  &
											 orientation, xc, yc                 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips real-valued array by elliptical obscuration centered at (xc,yc).
   ! majorDiam		Major axis of ellipse
   ! minorDiam		Minor axis
   ! orientation	Angle (counter-clockwise) between grid x-axis
   !					and major axis of ellipse
   !     17 Sept '02		Based on circular clip routine
   !     22 Nov  '02    Calls revised EllipticalIteration;
   !                    corrected calculation of the sign of edgeDist
   !     23 July '03     Either xc or yc (or both) may be omitted from input
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : R_2d

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr), DIMENSION(:),   INTENT(IN)     :: x, y
   REAL(pr),                 INTENT(IN)     :: majorDiam, minorDiam
   REAL(pr),                 INTENT(IN)     :: orientation
   REAL(pr),                 INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, PRESENT, MATMUL, SIGN

   ! Local variables
   REAL(pr), PARAMETER :: ZERO = 0.0_pr
   REAL(pr) :: dx, dy
   REAL(pr) :: dr
   REAL(pr) :: a, a_inner, a_outer
   REAL(pr) :: b, b_inner, b_outer
   REAL(pr) :: baSqr
   REAL(pr) :: xo, yo
   REAL(pr) :: xe, ye
   REAL(pr) :: edgeDist
   INTEGER  :: Ny, Nx, j, k
   REAL(pr), DIMENSION(2)   :: xy
   REAL(pr), DIMENSION(2,2) :: rotMatrix
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, minorDiam, xo, yo ) ) RETURN

   dr = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   a  = 0.5_pr * majorDiam
   b  = 0.5_pr * minorDiam
   a_inner = a - 2.0_pr*dr
   a_outer = a + 2.0_pr*dr
   b_inner = b - 2.0_pr*dr
   b_outer = b + 2.0_pr*dr

   baSqr = (b/a)**2

	rotMatrix = R_2d( -orientation )	! Transforms (x,y) -> (x',y'), a frame
												! in which x' and y' are along the major
												! and minor axes of the ellipse

   DO k = 1,Nx
   	DO j = 1,Ny
   		xy = MATMUL( rotMatrix, (/ x(k)-xo, y(j)-yo /) )

         IF ( (xy(1)/a_inner)**2 + (xy(2)/b_inner)**2 < 1.0_pr ) THEN
				! This point is completely inside the obscuration
            array(j,k) = ZERO
         ELSE IF ( (xy(1)/a_outer)**2 + (xy(2)/b_outer)**2 < 1.0_pr ) THEN
         	! This point is in the "fringe" region; do edge smoothing
         	CALL EllipticalIteration( xy(1), xy(2), a, b, xe, ye )
				edgeDist = SQRT( (xy(1)-xe)**2 + (xy(2)-ye)**2 )
				edgeDist = SIGN( edgeDist, 1.0_pr - ((xy(1)/a)**2 + (xy(2)/b)**2) )
				array(j,k) = SoftClip( array(j,k), edgeDist/dr )
         END IF

      END DO
   END DO

   END SUBROUTINE EllipticalObscReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE EllipticalObscCmplx( array, x, y, majorDiam, minorDiam,  &
											  orientation, xc, yc                 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips complex-valued array by elliptical obscuration centered at (xc,yc).
   ! majorDiam		Major axis of ellipse
   ! minorDiam		Minor axis
   ! orientation	Angle (counter-clockwise) between grid x-axis
   !					and major axis of ellipse
   !     17 Sept '02		Based on circular clip routine
   !     22 Nov  '02    Calls revised EllipticalIteration;
   !                    corrected calculation of the sign of edgeDist
   !     23 July '03    Either xc or yc (or both) may be omitted from input
   !     14 Oct  '03    ZERO declared COMPLEX (rather than REAL)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : R_2d

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x, y
   REAL(pr),                     INTENT(IN)     :: majorDiam, minorDiam
   REAL(pr),                     INTENT(IN)     :: orientation
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC SQRT, SIZE, PRESENT, MATMUL, SIGN

   ! Local variables
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )
   REAL(pr) :: dx, dy
   REAL(pr) :: dr
   REAL(pr) :: a, a_inner, a_outer
   REAL(pr) :: b, b_inner, b_outer
   REAL(pr) :: baSqr
   REAL(pr) :: xo, yo
   REAL(pr) :: xe, ye
   REAL(pr) :: edgeDist
   INTEGER  :: Ny, Nx, j, k
   REAL(pr), DIMENSION(2)   :: xy
   REAL(pr), DIMENSION(2,2) :: rotMatrix
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   IF ( PRESENT( yc ) ) yo = yc

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsideCircle( grid, minorDiam, xo, yo ) ) RETURN

   dr = SQRT( 0.5_pr * ( dx**2 + dy**2 ) ) * SMOOTHING_HALF_WIDTH
   a  = 0.5_pr * majorDiam
   b  = 0.5_pr * minorDiam
   a_inner = a - 2.0_pr*dr
   a_outer = a + 2.0_pr*dr
   b_inner = b - 2.0_pr*dr
   b_outer = b + 2.0_pr*dr

   baSqr = (b/a)**2

	rotMatrix = R_2d( -orientation )	! Transforms (x,y) -> (x',y'), a frame
												! in which x' and y' are along the major
												! and minor axes of the ellipse

   DO k = 1,Nx
   	DO j = 1,Ny
   		xy = MATMUL( rotMatrix, (/ x(k)-xo, y(j)-yo /) )

         IF ( (xy(1)/a_inner)**2 + (xy(2)/b_inner)**2 < 1.0_pr ) THEN
				! This point is completely inside the obscuration
            array(j,k) = ZERO
         ELSE IF ( (xy(1)/a_outer)**2 + (xy(2)/b_outer)**2 < 1.0_pr ) THEN
         	! This point is in the "fringe" region; do edge smoothing
         	CALL EllipticalIteration( xy(1), xy(2), a, b, xe, ye )
				edgeDist = SQRT( (xy(1)-xe)**2 + (xy(2)-ye)**2 )
				edgeDist = SIGN( edgeDist, 1.0_pr - ((xy(1)/a)**2 + (xy(2)/b)**2) )
				array(j,k) = SoftClip( array(j,k), edgeDist/dr )
         END IF

      END DO
   END DO

   END SUBROUTINE EllipticalObscCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE EllipticalIteration( x, y, a, b, xe, ye )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Finds point (xe,ye) on ellipse with normal that passes through (x,y).
   ! a = axis of ellipse along x
   ! b = axis of ellipse along y
	!
   ! Called by routines that clip by elliptical apertures or obscurations
   !     16 Sept '02
   !     19 Sept '02    Calls MakeArray
   !     22 Nov  '02    Revised input and improved initial guess for (xe,ye)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   REAL(pr), INTENT(IN)  :: x,  y
   REAL(pr), INTENT(IN)  :: a,  b
   REAL(pr), INTENT(OUT) :: xe, ye

   INTRINSIC RESHAPE, MATMUL, ABS, EPSILON

   ! Local variables
   REAL(pr) :: E, dEdx, dEdy
   REAL(pr) :: N, dNdx, dNdy
   REAL(pr) :: r
   REAL(pr) :: xi, eta
   REAL(pr) :: baSqr
   REAL(pr) :: xo, yo
   REAL(pr) :: dxo, dyo
   REAL(pr) :: det
   INTEGER  :: iteration
   REAL(pr), DIMENSION(2)   :: vector
   REAL(pr), DIMENSION(2,2) :: matrix
   INTEGER,  PARAMETER      :: ITER_MAX = 50
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   baSqr = (b/a)**2
   r     = SQRT( x**2 + y**2 )
   xi    = x/a
   eta   = y/b

	! First guess at solution is the radial projection of (x,y) onto ellipse:
	xo = x/r
	yo = y/r

	! Newton-Raphson iteration

	DO iteration = 1,ITER_MAX

		E = xo**2 + yo**2 - 1.0_pr				! Defines points on the ellipse.
		N = baSqr*xo*(eta-yo) - yo*(xi-xo)	! Defines normal to ellipse,
														! passing through (x,y)

		! Partial derivatives

		dEdx =  2.0_pr*xo
		dEdy =  2.0_pr*yo
		dNdx =  baSqr*(eta-yo) + yo
		dNdy = -baSqr*xo       - (xi-xo)

		! Solve simultaneous equations for corrections dxo and dyo

		matrix = MakeArray( 2, 2, (/  dNdy, -dEdy,   &
											  -dNdx,  dEdx /) )
		det    = dEdx*dNdy - dEdy*dNdx

		vector = (/ -E, -N /)
		vector = MATMUL( matrix, vector )/det

		dxo = vector(1)
		dyo = vector(2)
		xo = xo + dxo
		yo = yo + dyo

		IF ( ABS( dxo ) <= EPSILON( 1.0_pr )  .AND.  &
			  ABS( dyo ) <= EPSILON( 1.0_pr )         ) EXIT
	END DO

   IF ( iteration > ITER_MAX )  &
		CALL ErrorExit( "EllipticalIteration: Newton-Raphson did not converge" )

   xe = a*xo
   ye = b*yo

   END SUBROUTINE EllipticalIteration



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE PolygonClipReal( array, x, y, poly )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips a real-valued array by an aperture which is described by a convex
   ! polygon. For analytic geometry of convex polygons, see
   ! Angell & Griffith, "High Resolution Computer Graphics Using FORTRAN 77,"
   ! John Wiley & Sons (1987), pp.48-53.
   ! The edgeDist used here is positive outside the aperture, negative inside,
   ! and zero on the boundary.

   !     26 Oct  '00
   !      9 Mar  '01    "clockwise" and "notConvex" defined in Polygons module
   !     11 Mar  '01    Calls EdgeDistance function
   !      4 Dec  '01    Clear polyAp at end of calculation
   !      8 Feb  '02    Tests for aperture being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr), DIMENSION(:),   INTENT(IN)     :: x, y
   TYPE(polygon),            INTENT(IN)     :: poly

   INTRINSIC SIZE

   ! Local variables
   REAL(pr), PARAMETER :: ZERO = 0.0_pr
   INTEGER             :: orient
   TYPE(polygon)           :: polyAp
   REAL(pr)                :: dx, dy
   REAL(pr), DIMENSION(2)  :: newVertex, v, vp1
   REAL(pr)                :: wX, wY              ! Smoothing half-widths
   REAL(pr)                :: edgeDist
   INTEGER                 :: Nv, n, Nx, k, Ny, j
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nv = NumberOfVertices( poly )
   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsidePolygon( grid, poly ) ) RETURN

   polyAp = poly

   orient = PolygonOrientation( polyAp )
   IF ( orient == notConvex  .OR.  orient == clockwise ) &
      CALL ErrorExit( "PolygonClip: input polygon error" )

   ! Fill out the vertices of polyAp by repeating the first vertex
   ! to explicitly make a closed polygon

   newVertex = GetVertex( 1, poly )
   CALL AddVertex( polyAp, newVertex )

   wX = dx * SMOOTHING_HALF_WIDTH
   wY = dy * SMOOTHING_HALF_WIDTH

   vp1 = GetVertex( 1, polyAp )
   DO n = 1,Nv

      v  = vp1
      vp1 = GetVertex( n+1, polyAp )

      DO k = 1,Nx
      	DO j = 1,Ny

            edgeDist = EdgeDistance( x(k), y(j), v, vp1, wX, wY )
            IF (  edgeDist >= 1.0_pr ) THEN
               array(j,k) = ZERO
            ELSE IF ( edgeDist > -1.0_pr ) THEN
               array(j,k) = SoftClip( array(j,k), edgeDist )
            END IF

         END DO
      END DO

   END DO

   CALL ClearPolygon( polyAp )

   END SUBROUTINE PolygonClipReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE PolygonClipCmplx( array, x, y, poly )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips a complex-valued array by an aperture which is described by a convex
   ! polygon. For analytic geometry of convex polygons, see
   ! Angell & Griffith, "High Resolution Computer Graphics Using FORTRAN 77,"
   ! John Wiley & Sons (1987), pp.48-53.
   ! The edgeDist used here is positive outside the aperture, negative inside,
   ! and zero on the boundary.

   !     26 Oct  '00
   !      9 Mar  '01    "clockwise" and "notConvex" defined in Polygons module
   !     11 Mar  '01    Calls EdgeDistance function
   !      4 Dec  '01    Clear polyAp at end of calculation
   !      8 Feb  '02    Tests for aperture being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x,  y
   TYPE(polygon),                INTENT(IN)     :: poly

   INTRINSIC SIZE, REAL, AIMAG, CMPLX

   ! Local variables
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )
   INTEGER                 :: orient
   TYPE(polygon)           :: polyAp
   REAL(pr)                :: dx, dy
   REAL(pr), DIMENSION(2)  :: newVertex, v, vp1
   REAL(pr)                :: wX,    wY              ! Smoothing half-widths
   REAL(pr)                :: edgeDist
   INTEGER                 :: Nv, n, Nx, k, Ny, j
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nv = NumberOfVertices( poly )
   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsidePolygon( grid, poly ) ) RETURN

   polyAp = poly

   orient = PolygonOrientation( polyAp )
   IF ( orient == notConvex  .OR.  orient == clockwise ) &
      CALL ErrorExit( "PolygonClip: input polygon error" )

   ! Fill out the vertices of polyAp by repeating the first vertex
   ! to explicitly make a closed polygon

   newVertex = GetVertex( 1, poly )
   CALL AddVertex( polyAp, newVertex )

   wX = dx * SMOOTHING_HALF_WIDTH
   wY = dy * SMOOTHING_HALF_WIDTH

   vp1 = GetVertex( 1, polyAp )
   DO n = 1,Nv

      v  = vp1
      vp1 = GetVertex( n+1, polyAp )

      DO k = 1,Nx
      	DO j = 1,Ny

            edgeDist = EdgeDistance( x(k), y(j), v, vp1, wX, wY )
            IF (  edgeDist >= 1.0_pr ) THEN
               array(j,k) = ZERO
            ELSE IF ( edgeDist > -1.0_pr ) THEN
               array(j,k) = SoftClip( array(j,k), edgeDist )
            END IF

         END DO
      END DO

   END DO

   CALL ClearPolygon( polyAp )

   END SUBROUTINE PolygonClipCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE PolygonObscReal( array, x, y, poly )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips a real-valued array by an obscuration which is described by a convex
   ! polygon. For analytic geometry of convex polygons, see
   ! Angell & Griffith, "High Resolution Computer Graphics Using FORTRAN 77,"
   ! John Wiley & Sons (1987), pp.48-53.
   ! The edgeDist used here is positive outside the obscuration, negative inside,
   ! and zero on the boundary.

   !     26 Oct  '00
   !      9 Mar  '01    "clockwise" and "notConvex" defined in Polygons module
   !     11 Mar  '01    Test all edges for each grid point
   !      4 Dec  '01    Clear polyAp at end of calculation
   !      8 Feb  '02    Tests for obscuration being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr), DIMENSION(:),   INTENT(IN)     :: x,  y
   TYPE(polygon),            INTENT(IN)     :: poly

   INTRINSIC SIZE

   ! Local variables
   INTEGER                 :: orient
   TYPE(polygon)           :: polyAp
   REAL(pr)                :: dx, dy
   REAL(pr), DIMENSION(2)  :: newVertex, v, vp1
   REAL(pr)                :: factor
   REAL(pr)                :: wX, wY              ! Smoothing half-widths
   INTEGER                 :: Nv, n, Nx, k, Ny, j
   REAL(pr)                :: edgeDist
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nv = NumberOfVertices( poly )
   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsidePolygon( grid, poly ) )  &
      CALL ErrorExit( "PolygonObscReal: beam is completely obscured" )

   polyAp = poly

   orient = PolygonOrientation( polyAp )
   IF ( orient == notConvex  .OR.  orient == clockwise ) &
      CALL ErrorExit( "PolygonObsc: input polygon error" )

   ! Fill out the vertices of polyAp by repeating the first vertex
   ! to explicitly make a closed polygon

   newVertex = GetVertex( 1, poly )
   CALL AddVertex( polyAp, newVertex )

   wX = dx * SMOOTHING_HALF_WIDTH
   wY = dy * SMOOTHING_HALF_WIDTH

   ! Test each grid point against all polygon edges
   DO k = 1,Nx
   	DO j = 1,Ny

         factor = 1.0_pr
         ! For this grid point, calculate edge distances for all edges
         vp1 = GetVertex( 1, polyAp )
         DO n = 1,Nv
            v   = vp1
            vp1 = GetVertex( n+1, polyAp )
            edgeDist = EdgeDistance( x(k), y(j), v, vp1, wX, wY )
            IF( edgeDist > -1.0_pr )  factor = SoftClip( factor, edgeDist )
         END DO

         IF( factor > 0.0_pr ) array(j,k) = ( 1.0_pr - factor )*array(j,k)

      END DO
   END DO

   CALL ClearPolygon( polyAp )

   END SUBROUTINE PolygonObscReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE PolygonObscCmplx( array, x, y, poly )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips a wave by an obscuration which is described by a convex polygon.
   ! For analytic geometry of convex polygons, see
   ! Angell & Griffith, "High Resolution Computer Graphics Using FORTRAN 77,"
   ! John Wiley & Sons (1987), pp.48-53.
   ! The edgeDist used here is positive outside the obscuration, negative inside,
   ! and zero on the boundary.

   !     26 Oct  '00
   !      9 Mar  '01    "clockwise" and "notConvex" defined in Polygons module
   !     11 Mar  '01    Test all edges for each grid point
   !      4 Dec  '01    Clear polyAp at end of calculation
   !      8 Feb  '02    Tests for obscuration being larger than entire grid
   !		17 Sept '02		Calls SoftClip instead of TaperFunction
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x,  y
   TYPE(polygon),                INTENT(IN)     :: poly

   INTRINSIC SIZE, CMPLX, REAL, AIMAG

   ! Local variables
   INTEGER                 :: orient
   TYPE(polygon)           :: polyAp
   REAL(pr)                :: dx, dy
   REAL(pr), DIMENSION(2)  :: newVertex, v, vp1
   REAL(pr)                :: wX, wY              ! Smoothing half-widths
   REAL(pr)                :: factor, reA, imA
   INTEGER                 :: Nv, n, Nx, k, Ny, j
   REAL(pr)                :: edgeDist
   REAL(pr), DIMENSION(2,4) :: grid
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nv = NumberOfVertices( poly )
   Nx = SIZE( x )
   Ny = SIZE( y )
   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND = pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND = pr )

   grid = GridCornerArray( x(1), x(Nx), y(1), y(Ny) )
   IF ( GridIsInsidePolygon( grid, poly ) )  &
      CALL ErrorExit( "PolygonObscComplex: beam is completely obscured" )

   polyAp = poly

   orient = PolygonOrientation( polyAp )
   IF ( orient == notConvex  .OR.  orient == clockwise ) &
      CALL ErrorExit( "PolygonObsc: input polygon error" )

   ! Fill out the vertices of polyAp by repeating the first vertex
   ! to explicitly make a closed polygon

   newVertex = GetVertex( 1, poly )
   CALL AddVertex( polyAp, newVertex )

   wX = dx * SMOOTHING_HALF_WIDTH
   wY = dy * SMOOTHING_HALF_WIDTH

   ! Test each grid point against all polygon edges
   DO k = 1,Nx
   	DO j = 1,Ny

         factor = 1.0_pr
         ! For this grid point, calculate edge distances for all edges
         vp1 = GetVertex( 1, polyAp )
         DO n = 1,Nv
            v   = vp1
            vp1 = GetVertex( n+1, polyAp )
            edgeDist = EdgeDistance( x(k), y(j), v, vp1, wX, wY )
            IF( edgeDist > -1.0_pr )  factor = SoftClip( factor, edgeDist )
         END DO

         IF( factor > 0.0_pr ) THEN
            factor = 1.0_pr - factor
            reA    = factor * REAL( array(j,k) )
            imA    = factor * AIMAG( array(j,k) )
            array(j,k) = CMPLX( reA, imA, KIND=prc )
         END IF

      END DO
   END DO

   CALL ClearPolygon( polyAp )

   END SUBROUTINE PolygonObscCmplx


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION EdgeDistance( x, y, p1, p2, dx, dy ) RESULT( distance )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns the perpendicular distance from a point (x,y) to a directed line
   ! defined by point p1 = (x1,y1) and p2 = (x2,y2). The optional input dx and
   ! dy are scaling parameters for x and y coordinates, respectively. See
   ! Angell & Griffith, "High Resolution Computer Graphics Using FORTRAN 77,"
   ! John Wiley & Sons (1987), pp.48-53.

   !     11 Mar  '01    R.S. Benson
   !     29 May  '01    Revised
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),               INTENT(IN) :: x,  y
   REAL(pr), DIMENSION(2), INTENT(IN) :: p1, p2
   REAL(pr),     OPTIONAL, INTENT(IN) :: dx, dy

   REAL(pr) :: distance

   INTRINSIC PRESENT, SQRT

   ! Local variables
   REAL(pr) :: xDiff, yDiff, deltaX, deltaY
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xDiff = p2(1) - p1(1)
   yDiff = p2(2) - p1(2)

   deltaX = x - p1(1)
   deltaY = y - p1(2)

   IF( PRESENT( dx )  .AND.  PRESENT( dy ) ) THEN
      xDiff = xDiff/dx
      yDiff = yDiff/dy
      deltaX = deltaX/dx
      deltaY = deltaY/dy
   END IF

   distance = ( yDiff*deltaX - xDiff*deltaY )/SQRT( xDiff**2 + yDiff**2 )

   END FUNCTION EdgeDistance



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION GridCornerArray( xmin, xmax, ymin, ymax ) RESULT( gridCorners )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Creates a 2x4 array of coordinates of corners of a rectangular region

   ! xmin, xmax = extreme values of x-coordinate
   ! ymin, ymax = extreme values of y-coordinate
   !      8 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   REAL(pr), INTENT(IN)     :: xmin, xmax, ymin, ymax

   REAL(pr), DIMENSION(2,4) :: gridCorners

   ! Local variables
   REAL(pr), DIMENSION(8) :: vector
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   vector = (/ xmin, xmax, xmax, xmin, &
               ymin, ymin, ymax, ymax /)
   gridCorners = MakeArray( 2, 4, vector )

   END FUNCTION GridCornerArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION GridIsInsideCircle( gridCorners, diam, xc, yc ) RESULT( answer )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Tests four corners of rectangular grid; returns answer = .TRUE. if
   ! they are all inside the circle

   ! gridCorners = array of coordinates of corners of rectangle
   ! diam        = size of circle
   ! xc, yc      = coordinates of center of circle

   !      8 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(2,4), INTENT(IN) :: gridCorners
   REAL(pr),                 INTENT(IN) :: diam, xc, yc

   LOGICAL :: answer

   ! Local variables
   REAL(pr) :: x, y
   REAL(pr) :: radSqr
   INTEGER  :: corner
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   radSqr = (0.5_pr*diam)**2

   answer = .TRUE.
   DO corner = 1,4
      x = gridCorners(1,corner)
      y = gridCorners(2,corner)
      IF ( ( (x-xc)**2 + (y-yc)**2 - radSqr ) >= 0.0_pr ) answer = .FALSE.
   END DO

   END FUNCTION GridIsInsideCircle



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION GridIsInsidePolygon( gridCorners, poly ) RESULT( answer )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Tests four corners of rectangular grid; returns answer = .TRUE. if
   ! they are all inside the polygon

   ! gridCorners = array of coordinates of corners of rectangle
   ! poly        = a convex polygon

   !      8 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Polygons

   IMPLICIT NONE

   REAL(pr), DIMENSION(2,4), INTENT(IN) :: gridCorners
   TYPE(polygon),            INTENT(IN) :: poly

   LOGICAL :: answer

   ! Local variables
   REAL(pr), DIMENSION(2) :: p1, p2    ! Points defining one side of polygon
   REAL(pr) :: x,  y
   INTEGER  :: Nv, n
   INTEGER  :: corner
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nv = NumberOfVertices( poly )

   answer = .TRUE.
   p1 = GetVertex( Nv, poly )
   DO n = 1,Nv
      p2 = GetVertex( n, poly )
      DO corner = 1,4
         x = gridCorners(1,corner)
         y = gridCorners(2,corner)
         IF ( EdgeDistance( x, y, p1, p2 ) >= 0.0_pr ) answer = .FALSE.
      END DO
      p1 = p2
   END DO

   END FUNCTION GridIsInsidePolygon



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION SoftClipReal( inputValue, r ) RESULT( outputValue )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips using a function that varies smoothly and monotonically between
   ! 1 and 0, as r varies from -1 to + 1. The particular function
   ! used here is the convolution of a step function with
   ! a triangle function.

   ! inputValue = Real-valued sample value to be clipped
   ! r = distance from aperture edge, in units of the smoothing half-width
   !     17 Sept '02		Based on TaperFunction
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN) :: inputValue
   REAL(pr), INTENT(IN) :: r
   REAL(pr)             :: outputValue

   INTRINSIC ABS

   ! Local variables
   REAL(pr) :: x
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   x = 0.5_pr * ( 1.0_pr - ABS( r ) )**2

   outputValue = inputValue

   IF ( r >= 1.0_pr ) THEN
      outputValue = 0.0_pr
   ELSE IF ( r >= 0.0_pr ) THEN
      outputValue = x * inputValue
   ELSE IF ( r > -1.0_pr ) THEN
      outputValue = (1.0_pr - x) * inputValue
   END IF

   END FUNCTION SoftClipReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION SoftClipCmplx( inputValue, r ) RESULT( outputValue )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Clips using a function that varies smoothly and monotonically between
   ! 1 and 0, as r varies from -1 to + 1. The particular function
   ! used here is the convolution of a step function with
   ! a triangle function.

   ! inputValue = Complex-valued sample value to be clipped
   ! r = distance from aperture edge, in units of the smoothing half-width
   !     17 Sept '02		Based on TaperFunction
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), INTENT(IN) :: inputValue
   REAL(pr),     INTENT(IN) :: r
   COMPLEX(prc)             :: outputValue

   INTRINSIC ABS, REAL, AIMAG, CMPLX

   ! Local variables
   REAL(pr) :: x
   REAL(pr) :: realPart, imagPart
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   x = 0.5_pr * ( 1.0_pr - ABS( r ) )**2

   realPart = REAL( inputValue )
   imagPart = AIMAG( inputValue )

   IF ( r >= 1.0_pr ) THEN
      realPart = 0.0_pr
      imagPart = 0.0_pr
   ELSE IF ( r >= 0.0_pr ) THEN
      realPart = x * realPart
      imagPart = x * imagPart
   ELSE IF ( r > -1.0_pr ) THEN
      realPart = (1.0_pr - x) * realPart
      imagPart = (1.0_pr - x) * imagPart
   END IF

	outputValue = CMPLX( realPart, imagPart, KIND=prc )

   END FUNCTION SoftClipCmplx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! ADDITIONAL SUBROUTINES
!!!!! David Marx
!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! Attempt at making a mask using a window for low diffraction sidelobes
!! April 16, 2004
!! 
SUBROUTINE WindowClip(beamamp, x, y, dx, dy, winT, winA, xc, yc, rotAngle, obsc)

IMPLICIT NONE
   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: beamamp
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x, y
   REAL(pr),                     INTENT(IN)     :: dx, dy 
   REAL(pr),                     INTENT(IN)     :: winT, winA
   REAL(pr),                     INTENT(IN)     :: xc, yc, rotAngle
   LOGICAL,                      INTENT(IN)     :: obsc


   ! Local Variables
   INTEGER                 :: Nx, Ny, i, j
   REAL(pr)                :: xp, yp
   REAL(pr)                :: redge, rnorm, rtmp
   REAL(pr)                :: wx, dwx
   REAL(pr)                :: obscsign = 1.0_pr
   REAL(pr), PARAMETER     :: pi = 3.14159265358979_pr
   COMPLEX(prc), PARAMETER :: ZERO = ( 0.0_pr, 0.0_pr )

   !!!!!!!!!!!!!!!!!!!!!!!!
   Nx = SIZE(x)
   Ny = SIZE(y)

   IF ( obsc ) obscsign = -1.0_pr

   rnorm = SQRT(dx*dx + dy*dy)

   DO i = 1,Nx

      Do j = 1,Ny

         !! first apply coordinate transformation
         xp = (x(i)-xc)*cos(rotAngle) + (y(j)-yc)*sin(rotAngle)
         yp =-(x(i)-xc)*sin(rotAngle) + (y(j)-yc)*cos(rotAngle)

         ! test if x is within the window
         IF ( xp > -winT/2 .AND. xp < winT/2 ) THEN

            ! the window function and its derivative (slope) at x(i)
             wx = 0.42_pr + 0.5_pr*COS(2.0_pr*pi*xp/winT) + 0.08_pr*COS(4.0_pr*pi*xp/winT)
             wx = winA * wx
             dwx = -(pi/winT)*SIN(2.0_pr*pi*xp/winT) - 0.32_pr*(pi/winT)*SIN(4.0_pr*pi*xp/winT)
             dwx = winA*dwx

             ! edge distance at this point of the window
             rtmp = obscsign * (1/SQRT(1+dwx*dwx))

             redge = ( ABS(yp)-wx )*rtmp/rnorm

             IF ( redge >= 1.0_pr ) THEN
                beamamp(j,i) = ZERO
             ELSE IF ( redge > -1.0_pr ) THEN
                beamamp(j,i) = SoftClip( beamamp(j,i), redge )
             END IF ! otherwise x,y is inside window and beamamp is unchanged
                

 

          ELSE  ! x is before or after the window

             IF ( .NOT.obsc ) beamamp(j,i) = ZERO

          END IF ! if x is within the window

      END DO ! j = 1,Ny

   END DO ! i = 1,Nx
         

END SUBROUTINE WindowClip




END MODULE Clipping_Routines
