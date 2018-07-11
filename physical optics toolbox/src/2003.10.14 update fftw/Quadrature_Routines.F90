!     Last change:  RSB   1 Mar 2002    3:11 pm
MODULE Quadrature_Routines

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

   INTRINSIC RESHAPE

   PUBLIC
   PRIVATE :: Quad2dReal,            Quad2dCmplx,         &
              Quad2dRectReal,        Quad2dRectCmplx,     &
              Quad2dCircReal,        Quad2dCircCmplx,     &
              Quad2dFullgridReal,    Quad2dFullgridCmplx, &
              Quad1DReal,            Quad1DCmplx,         &
              IntegrateRealOverArea, IntegrateCmpxOverArea

   INTERFACE IntegrateOverArea
      MODULE PROCEDURE IntegrateRealOverArea, &
                       IntegrateCmpxOverArea
   END INTERFACE

   INTERFACE Quad2d
      MODULE PROCEDURE Quad2dReal, &
                       Quad2dCmplx
   END INTERFACE

   INTERFACE Quad1D
      MODULE PROCEDURE Quad1DReal, &
                       Quad1DCmplx
   END INTERFACE

   ! Matrix used by Quad1DReal and Quad1DCmplx
   INTEGER,  DIMENSION(2),   PARAMETER, PRIVATE :: rowOrder = (/ 2, 1 /)
   REAL(pr), DIMENSION(4,4), PARAMETER, PRIVATE ::                             &
   	MFIT = (1.0_pr/24.0_pr)*RESHAPE( (/ 0.0_pr, 24.0_pr,  0.0_pr,  0.0_pr,   &
      	                                 2.0_pr, -6.0_pr,  6.0_pr, -2.0_pr,   &
      	                                 4.0_pr, -8.0_pr,  4.0_pr,  0.0_pr,   &
      	                                -1.0_pr,  3.0_pr, -3.0_pr,  1.0_pr /),&
      	                                 (/4,4/), ORDER=rowOrder )

   ! Constants used by CurveFitCoeffsReal and CurveFitCoeffsCmplx
   REAL(pr), PARAMETER, PRIVATE :: c23 = 2.0_pr/3.0_pr
   REAL(pr), PARAMETER, PRIVATE :: c25 = 2.0_pr/5.0_pr
   REAL(pr), PARAMETER, PRIVATE :: c27 = 2.0_pr/7.0_pr

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION IntegrateRealOverArea( array, dx, dy,                 &
                                   xc, yc, shape, xSize, ySize )  &
                                 RESULT( integral )

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Does 2-d integration of input real-valued array over specified area. The
   ! optional input "ySize" is needed only for shape = "rectangle." Integration
   ! is a simple trapezoid approximation.

   !     20 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : CoordinateVector
   USE Clipping_Routines
   USE Polygons

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN) :: array
   REAL(pr),                 INTENT(IN) :: dx, dy
   REAL(pr),                 INTENT(IN) :: xc, yc
   CHARACTER( LEN=* ),       INTENT(IN) :: shape
   REAL(pr),                 INTENT(IN) :: xSize
   REAL(pr),                 INTENT(IN), OPTIONAL :: ySize

   REAL(pr) :: integral

   INTRINSIC SIZE, PRESENT, SUM

   ! Local variables
   REAL(pr), DIMENSION(SIZE( array, 1),SIZE( array,2 )) :: TwoDarray
   REAL(pr), DIMENSION(SIZE( array, 1)) :: y
   REAL(pr), DIMENSION(SIZE( array, 2)) :: x
   REAL(pr)      :: width, angle
   TYPE(polygon) :: area
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   x = CoordinateVector( SIZE( array, 2 ), dx )
   y = CoordinateVector( SIZE( array, 1 ), dy )
   TwoDarray = array

   SELECT CASE ( shape )
      CASE ( "circle" )
         CALL CircularClip( TwoDarray, x, y, xSize, xc, yc )
      CASE ( "square", "rectangle" )
         IF ( shape == "square" ) THEN
            width = xSize
         ELSE IF ( shape == "rectangle"  .AND.  PRESENT( ySize ) ) THEN
            width = ySize
         ELSE
            CALL ErrorExit( "IntegrateRealOverArea: missing rectangle data" )
         END IF
         angle = 0.0_pr
         area  = MakeRectangle( xSize, width, angle, xc, yc )
         CALL PolygonClip( TwoDarray, x, y, area )
      CASE DEFAULT
         CALL ErrorExit( "IntegrateRealOverArea: shape error" )
   END SELECT

   integral = (dx*dy)*SUM( TwoDarray )

   END FUNCTION IntegrateRealOverArea



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION IntegrateCmpxOverArea( array, dx, dy,                 &
                                   xc, yc, shape, xSize, ySize )  &
                                 RESULT( integral )

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Does 2-d integration of input complex-valued array over specified area. The
   ! optional input "ySize" is needed only for shape = "rectangle." Integration
   ! is a simple trapezoid approximation.

   !     20 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines, ONLY : CoordinateVector
   USE Clipping_Routines
   USE Polygons

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: array
   REAL(pr),                     INTENT(IN) :: dx, dy
   REAL(pr),                     INTENT(IN) :: xc, yc
   CHARACTER( LEN=* ),           INTENT(IN) :: shape
   REAL(pr),                     INTENT(IN) :: xSize
   REAL(pr),                     INTENT(IN), OPTIONAL :: ySize

   COMPLEX(prc) ::   integral

   INTRINSIC SIZE, PRESENT, SUM

   ! Local variables
   COMPLEX(prc), DIMENSION(SIZE( array, 1),SIZE( array,2 )) :: TwoDarray
   REAL(pr),     DIMENSION(SIZE( array, 1)) :: y
   REAL(pr),     DIMENSION(SIZE( array, 2)) :: x
   REAL(pr)      :: width, angle
   TYPE(polygon) :: area
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   x = CoordinateVector( SIZE( array, 2 ), dx )
   y = CoordinateVector( SIZE( array, 1 ), dy )
   TwoDarray = array

   SELECT CASE ( shape )
      CASE ( "circle" )
         CALL CircularClip( TwoDarray, x, y, xSize, xc, yc )
      CASE ( "square", "rectangle" )
         IF ( shape == "square" ) THEN
            width = xSize
         ELSE IF ( shape == "rectangle"  .AND.  PRESENT( ySize ) ) THEN
            width = ySize
         ELSE
            CALL ErrorExit( "IntegrateCmpxOverArea: missing rectangle data" )
         END IF
         angle = 0.0_pr
         area  = MakeRectangle( xSize, width, angle, xc, yc )
         CALL PolygonClip( TwoDarray, x, y, area )
      CASE DEFAULT
         CALL ErrorExit( "IntegrateCmpxOverArea: shape error" )
   END SELECT

   integral = (dx*dy)*SUM( TwoDarray )

   END FUNCTION IntegrateCmpxOverArea



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dReal( array, x, y, shape, size1, size2, xc, yc )  &
                      RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates a two-dimensional real-valued quadrature over a specified area.
   ! The shape may be circle, square, rectangle, or "infinite" (entire grid);
   ! if the shape parameter is not present, "infinite" is assumed

   !      28 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN) :: array
   REAL(pr), DIMENSION(:),   INTENT(IN) :: x, y
   CHARACTER(LEN=*),         INTENT(IN), OPTIONAL :: shape
   REAL(pr),                 INTENT(IN), OPTIONAL :: size1, size2
   REAL(pr),                 INTENT(IN), OPTIONAL :: xc, yc

   REAL(pr) :: quad

   INTRINSIC PRESENT

   ! Local variables
   CHARACTER(LEN=9) :: areaShape
   REAL(pr)         :: xo, yo
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   areaShape = "infinite"
   IF ( PRESENT( shape ) ) areaShape = shape
   xo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   yo = 0.0_pr
   IF ( PRESENT( yc ) ) yo = yc

   SELECT CASE ( areaShape )
      CASE( "circle", "Circle", "CIRCLE" )
         IF ( .NOT.PRESENT( size1 ) )  &
         	CALL ErrorExit( "Quad2d: diameter not specified" )
         quad = Quad2dCircReal( array, x, y, size1, xo, yo )
      CASE( "square", "Square", "SQUARE" )
         IF ( .NOT.PRESENT( size1 ) )  &
         	CALL ErrorExit( "Quad2d: size of square not specified" )
         quad = Quad2dRectReal( array, x, y, size1, size1, xo, yo )
      CASE( "rectangle", "Rectangle", "RECTANGLE" )
         IF ( .NOT.PRESENT( size1 )  .OR.  .NOT.PRESENT( size2 ) )  &
         	CALL ErrorExit( "Quad2d: rectangle size not specified" )
         quad = Quad2dRectReal( array, x, y, size1, size2, xo, yo )
      CASE( "infinite", "Infinite", "INFINITE" )
         quad = Quad2dFullgridReal( array, x, y )
      CASE DEFAULT
         CALL ErrorExit( "Quad2d: error specifying detector shape" )
   END SELECT

   END FUNCTION Quad2dReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dCmplx( array, x, y, shape, size1, size2, xc, yc )  &
                       RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates a two-dimensional complex-valued quadrature over a specified
   ! area. The shape may be circle, square, rectangle, or "infinite";
   ! if the shape parameter is not present, "infinite" (entire grid) is assumed.
   ! Vectors x and y are coordinate values for function values of the array.
   ! The center of the area is located at coordinates (xc,yc).

   !      28 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: array
   REAL(pr),     DIMENSION(:),   INTENT(IN) :: x, y
   CHARACTER(LEN=*),             INTENT(IN), OPTIONAL :: shape
   REAL(pr),                     INTENT(IN), OPTIONAL :: size1, size2
   REAL(pr),                     INTENT(IN), OPTIONAL :: xc, yc

   COMPLEX(prc) :: quad

   INTRINSIC PRESENT

   ! Local variables
   CHARACTER(LEN=9) :: areaShape
   REAL(pr)         :: xo, yo
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   areaShape = "infinite"
   IF ( PRESENT( shape ) ) areaShape = shape
   xo = 0.0_pr
   IF ( PRESENT( xc ) ) xo = xc
   yo = 0.0_pr
   IF ( PRESENT( yc ) ) yo = yc

   SELECT CASE ( areaShape )
      CASE( "circle", "Circle", "CIRCLE" )
         IF ( .NOT.PRESENT( size1 ) )  &
         	CALL ErrorExit( "Quad2d: diameter not specified" )
         quad = Quad2dCircCmplx( array, x, y, size1, xo, yo )
      CASE( "square", "Square", "SQUARE" )
         IF ( .NOT.PRESENT( size1 ) )  &
         	CALL ErrorExit( "Quad2d: size of square not specified" )
         quad = Quad2dRectCmplx( array, x, y, size1, size1, xo, yo )
      CASE( "rectangle", "Rectangle", "RECTANGLE" )
         IF ( .NOT.PRESENT( size1 )  .OR.  .NOT.PRESENT( size2 ) )  &
         	CALL ErrorExit( "Quad2d: rectangle size not specified" )
         quad = Quad2dRectCmplx( array, x, y, size1, size2, xo, yo )
      CASE( "infinite", "Infinite", "INFINITE" )
         quad = Quad2dFullgridCmplx( array, x, y )
      CASE DEFAULT
         CALL ErrorExit( "Quad2d: error specifying detector shape" )
   END SELECT

   END FUNCTION Quad2dCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dRectReal( TwoDfunc, x, y, xSize, ySize, xc, yc )  &
                              RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the integral of the input real-valued function over a
   ! rectangular area.

   !      2 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN) :: TwoDfunc
   REAL(pr), DIMENSION(:),   INTENT(IN) :: x, y
   REAL(pr),                 INTENT(IN) :: xSize, ySize     ! Rectangle size
   REAL(pr),                 INTENT(IN) :: xc, yc           ! Rectangle center

   REAL(pr) :: quad

   ! Local variables
   REAL(pr) :: xLower, xUpper, yLower, yUpper
   REAL(pr), DIMENSION(:), ALLOCATABLE :: yQuad
   INTEGER  :: k, k1, k2, j1, j2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xLower = xc - 0.5_pr*xSize
   xUpper = xc + 0.5_pr*xSize
   yLower = yc - 0.5_pr*ySize
   yUpper = yc + 0.5_pr*ySize

   k1 = FindLowerIndex( xLower, x )
   j1 = FindLowerIndex( yLower, y )
   k2 = FindUpperIndex( xUpper, x )
   j2 = FindUpperIndex( yUpper, y )

   ALLOCATE( yQuad(k1-2:k2+2) )

   ! First integrate along y, then along x

   DO k = k1-2,k2+2
      yQuad(k) = Quad1D( TwoDfunc(:,k), y, yLower, yUpper, j1, j2 )
   END DO

   quad = Quad1D( yQuad, x(k1-2:k2+2), xLower, xUpper )

   DEALLOCATE( yQuad)

   END FUNCTION Quad2dRectReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dRectCmplx( TwoDfunc, x, y, xSize, ySize, xc, yc )  &
                               RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the integral of the input complex-valued function over a
   ! rectangular area.

   !      2 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: TwoDfunc
   REAL(pr),     DIMENSION(:),   INTENT(IN) :: x, y
   REAL(pr),                     INTENT(IN) :: xSize, ySize ! Rectangle size
   REAL(pr),                     INTENT(IN) :: xc, yc       ! Rectangle center

   COMPLEX(prc) :: quad

   ! Local variables
   REAL(pr) :: xLower, xUpper, yLower, yUpper
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE ::   yQuad
   INTEGER  :: k, k1, k2, j1, j2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xLower = xc - 0.5_pr*xSize
   xUpper = xc + 0.5_pr*xSize
   yLower = yc - 0.5_pr*ySize
   yUpper = yc + 0.5_pr*ySize

   k1 = FindLowerIndex( xLower, x )
   j1 = FindLowerIndex( yLower, y )
   k2 = FindUpperIndex( xUpper, x )
   j2 = FindUpperIndex( yUpper, y )

   ALLOCATE( yQuad(k1-2:k2+2) )

   ! First integrate along y, then along x

   DO k = k1-2,k2+2
      yQuad(k) = Quad1D( TwoDfunc(:,k), y, yLower, yUpper, j1, j2 )
   END DO

   quad = Quad1D( yQuad, x(k1-2:k2+2), xLower, xUpper )

   DEALLOCATE( yQuad )

   END FUNCTION Quad2dRectCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dCircReal( TwoDfunc, x, y, diam, xc, yc ) RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the integral of the input real-valued function over a
   ! circular area.

   !      2 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN) :: TwoDfunc
   REAL(pr), DIMENSION(:),   INTENT(IN) :: x, y          ! Coordinate vectors
   REAL(pr),                 INTENT(IN) :: diam
   REAL(pr),                 INTENT(IN) :: xc, yc        ! Circle's center

   REAL(pr) :: quad

   INTRINSIC SQRT

   ! Local variables
   REAL(pr) :: radius
   REAL(pr) :: xLower, xUpper, yLower, yUpper
   INTEGER  :: k, k1, k2, j1, j2
   REAL(pr), DIMENSION(:), ALLOCATABLE :: yQuad
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   radius = 0.5_pr*diam

   xLower = xc - radius
   xUpper = xc + radius

   k1 = FindLowerIndex( xLower, x )
   k2 = FindUpperIndex( xUpper, x )

   ALLOCATE( yQuad(k1:k2) )

   ! First integrate along y, then along x

   DO k = k1,k2
      yLower = yc - SQRT( radius**2 - (x(k)-xc)**2 )
      yUpper = yc + SQRT( radius**2 - (x(k)-xc)**2 )
      yQuad(k) = Quad1D( TwoDfunc(:,k), y, yLower, yUpper )
   END DO

   quad = Xquad1dReal( yQuad, x(k1:k2), xLower, xUpper )

   DEALLOCATE( yQuad )

   CONTAINS

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FUNCTION Xquad1dReal( xFunc, x, xL, xU ) RESULT( quad )
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      !   3 Nov  '00
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      USE Math_Routines

      IMPLICIT NONE

      REAL(pr), DIMENSION(:), INTENT(IN) :: xFunc     ! Function values
      REAL(pr), DIMENSION(:), INTENT(IN) :: x         ! Independent variable
      REAL(pr),               INTENT(IN) :: xL, xU    ! Lower and upper limits

      REAL(pr) :: quad

      INTRINSIC SIZE, MODULO, SUM, SQRT

      ! Local variables
      REAL(pr) :: dx
      REAL(pr) :: qSimp, lowerEnd, upperEnd
      REAL(pr) :: y
      INTEGER  :: N, nStart, kStart
      REAL(pr), DIMENSION(3)          :: z, f, coeffs
      REAL(pr),             PARAMETER ::  CRIT = 1.0d-2
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      N  = SIZE( xFunc )
      dx = ( x(N) - x(1) )/REAL( N-1, KIND=pr )

      IF ( N < 4 ) CALL ErrorExit( "Xquad1dReal: too few points for quadrature" )

      nStart = 1
      IF ( MODULO( N-1, 2 ) /= 0 ) nStart = 2

      ! Do Simpson quadrature for interior points (nStart to N)

      qSimp = dx * ( xFunc(nStart) + xFunc(N) +                    &
                     4.0_pr*SUM( xFunc(nStart+1:N-1:2) ) +         &
                     2.0_pr*SUM( xFunc(nStart+2:N-2:2) )  ) / 3.0_pr

      ! Add contributions for ends of integration interval

      kStart = 1
      IF ( x(1)-xL < CRIT*dx ) kStart = 2
      z = x(kStart:kStart+2) - xL
      f = xFunc(kStart:kStart+2) / SQRT( z )
      coeffs = CalcCurveFitCoeffsReal( z, dx, f )
      y = x(nStart) - xL
      lowerEnd = y * SQRT( y ) * EvalPoly( y, coeffs )

      kStart = N
      IF ( xU-x(kStart) < CRIT*dx ) kStart = N-1
      z = xU - x(kStart:kStart-2:-1)
      f = xFunc(kStart:kStart-2:-1) / SQRT( z )
      coeffs = CalcCurveFitCoeffsReal( z, dx, f )
      y = xU - x(N)
      upperEnd = y * SQRT( y ) * EvalPoly( y, coeffs )

      quad = qSimp + lowerEnd + upperEnd

      END FUNCTION Xquad1dReal

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FUNCTION CalcCurveFitCoeffsReal( z, dz, f ) RESULT( coeffs )
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Returns coefficients for a curve fit of the form
      !           f(x) = sqrt(x) * ( a + bx + cx**2 )
      !
      !                  coeff(1) = (2/3)a
      !                  coeff(2) = (2/5)b
      !                  coeff(3) = (2/7)c

      !   3 Nov  '00
      !   1 Jun  '01    2/3, 2/5, and 2/7 are private module PARAMETERS
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      REAL(pr), DIMENSION(3), INTENT(IN) :: z, f
      REAL(pr),               INTENT(IN) :: dz

      REAL(pr), DIMENSION(3) :: coeffs

      ! Local variables (none)
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      coeffs(3) = ( f(3) - f(1) )/(2.0_pr*dz**2)
      coeffs(2) = ( f(3) - f(2) )/dz - ( z(3) + z(2) )*coeffs(3)
      coeffs(1) = f(1) - z(1)*( coeffs(2) + z(1)*coeffs(3) )

      coeffs(1) = c23*coeffs(1)
      coeffs(2) = c25*coeffs(2)
      coeffs(3) = c27*coeffs(3)

      END FUNCTION CalcCurveFitCoeffsReal

   END FUNCTION Quad2dCircReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dCircCmplx( TwoDfunc, x, y, diam, xc, yc ) RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the integral of the input complex-valued function over a
   ! circular area.

   !      2 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: TwoDfunc
   REAL(pr),     DIMENSION(:),   INTENT(IN) :: x, y         ! Coordinate vectors
   REAL(pr),                     INTENT(IN) :: diam
   REAL(pr),                     INTENT(IN) :: xc, yc       ! Circle's center

   COMPLEX(prc) :: quad

   INTRINSIC SQRT

   ! Local variables
   REAL(pr) :: radius
   REAL(pr) :: xLower, xUpper, yLower, yUpper
   INTEGER  :: k, k1, k2, j1, j2
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE ::   yQuad
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   radius = 0.5_pr*diam

   xLower = xc - radius
   xUpper = xc + radius

   k1 = FindLowerIndex( xLower, x )
   k2 = FindUpperIndex( xUpper, x )

   ALLOCATE( yQuad(k1:k2) )

   ! First integrate along y, then along x

   DO k = k1,k2
      yLower = yc - SQRT( radius**2 - (x(k)-xc)**2 )
      yUpper = yc + SQRT( radius**2 - (x(k)-xc)**2 )
      yQuad(k) = Quad1D( TwoDfunc(:,k), y, yLower, yUpper )
   END DO

   quad = Xquad1dCmplx( yQuad, x(k1:k2), xLower, xUpper )

   DEALLOCATE( yQuad )

   CONTAINS

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FUNCTION Xquad1dCmplx( xFunc, x, xL, xU ) RESULT( quad )
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      !   3 Nov  '00
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      USE Math_Routines

      IMPLICIT NONE

      COMPLEX(prc), DIMENSION(:), INTENT(IN) :: xFunc    ! Function values
      REAL(pr),     DIMENSION(:), INTENT(IN) :: x        ! Independent variable
      REAL(pr),                   INTENT(IN) :: xL, xU   ! Lower and upper limits

      COMPLEX(prc) :: quad

      INTRINSIC SIZE, MODULO, SUM, SQRT

      ! Local variables
      REAL(pr)     ::   dx
      COMPLEX(prc) ::   qSimp, lowerEnd, upperEnd
      REAL(pr)     :: y
      INTEGER      :: N, nStart, kStart
      REAL(pr),     DIMENSION(3)          :: z
      COMPLEX(prc), DIMENSION(3)          :: f, coeffs
      REAL(pr),                 PARAMETER :: CRIT = 1.0d-2
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      N  = SIZE( xFunc )
      dx = ( x(N) - x(1) )/REAL( N-1, KIND=pr )

      IF ( N < 4 ) CALL ErrorExit( "Xquad1dCmplx: too few points for quadrature" )

      nStart = 1
      IF ( MODULO( N-1, 2 ) /= 0 ) nStart = 2

      ! Do Simpson quadrature for interior points (nStart to N)

      qSimp = dx * ( xFunc(nStart) + xFunc(N) +                    &
                     4.0_pr*SUM( xFunc(nStart+1:N-1:2) ) +         &
                     2.0_pr*SUM( xFunc(nStart+2:N-2:2) )  ) / 3.0_pr

      ! Add contributions for ends of integration interval

      kStart = 1
      IF ( x(1)-xL < CRIT*dx ) kStart = 2
      z = x(kStart:kStart+2) - xL
      f = xFunc(kStart:kStart+2) / SQRT( z )
      coeffs = CalcCurveFitCoeffsCmplx( z, dx, f )
      y = x(nStart) - xL
      lowerEnd = y * SQRT( y ) * EvalPoly( y, coeffs )

      kStart = N
      IF ( xU-x(kStart) < CRIT*dx ) kStart = N-1
      z = xU - x(kStart:kStart-2:-1)
      f = xFunc(kStart:kStart-2:-1) / SQRT( z )
      coeffs = CalcCurveFitCoeffsCmplx( z, dx, f )
      y = xU - x(N)
      upperEnd = y * SQRT( y ) * EvalPoly( y, coeffs )

      quad = qSimp + lowerEnd + upperEnd

      END FUNCTION Xquad1dCmplx

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FUNCTION CalcCurveFitCoeffsCmplx( z, dz, f ) RESULT( coeffs )
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Returns coefficients for a curve fit of the form
      !           f(x) = sqrt(x) * ( a + bx + cx**2 )
      !
      !                  coeff(1) = (2/3)a
      !                  coeff(2) = (2/5)b
      !                  coeff(3) = (2/7)c

      !   3 Nov  '00
      !   1 Jun  '01    2/3, 2/5, and 2/7 are private module PARAMETERS
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      REAL(pr),     DIMENSION(3), INTENT(IN) :: z
      REAL(pr),                   INTENT(IN) :: dz
      COMPLEX(prc), DIMENSION(3), INTENT(IN) :: f

      COMPLEX(prc), DIMENSION(3) :: coeffs

      ! Local variables (none)
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      coeffs(3) = ( f(3) - f(1) )/(2.0_pr*dz**2)
      coeffs(2) = ( f(3) - f(2) )/dz - ( z(3) + z(2) )*coeffs(3)
      coeffs(1) = f(1) - z(1)*( coeffs(2) + z(1)*coeffs(3) )

      coeffs(1) = c23*coeffs(1)
      coeffs(2) = c25*coeffs(2)
      coeffs(3) = c27*coeffs(3)

      END FUNCTION CalcCurveFitCoeffsCmplx

   END FUNCTION Quad2dCircCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dFullgridReal( TwoDfunc, x, y ) RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the integral of the input real-valued function over the full
   ! two-dimensional grid, using successive 1-d quadratures. The function
   ! values are assumed to be uniformly-spaced.

   !      9 July '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN) :: TwoDfunc
   REAL(pr), DIMENSION(:),   INTENT(IN) :: x, y

   REAL(pr) :: quad

   INTRINSIC SIZE, REAL

   ! Local variables
   REAL(pr), DIMENSION(SIZE( TwoDfunc, 1 )) :: yFunc
   INTEGER  :: Nx, Ny, j
   REAL(pr) :: dx, dy
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( TwoDfunc, 2 )
   Ny = SIZE( TwoDfunc, 1 )

   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND=pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND=pr )

   ! First integrate along x, for each y-value...

   DO j = 1, Ny
      yFunc(j) = OneDquadReal( TwoDfunc(j,:), Nx )
   END DO

   ! ... then integrate along y

   quad = (dx*dy)*OneDquadReal( yFunc, Ny )

   CONTAINS

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FUNCTION OneDquadReal( func, N ) RESULT( quadValue )
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Uses a closed-interval, extended formula of order (1/N^4), from
      ! "Numerical Recipes in Fortran," Second Edition, equation (4.1.14).
      !   9 July '01
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      REAL(pr), DIMENSION(:), INTENT(IN) :: func
      INTEGER,                INTENT(IN) :: N

      REAL(pr)                           :: quadValue

      INTRINSIC SUM

      ! Local variables:
      REAL(pr), PARAMETER :: c1 = 3.0_pr/8.0_pr
      REAL(pr), PARAMETER :: c2 = 7.0_pr/6.0_pr
      REAL(pr), PARAMETER :: c3 = 23.0_pr/24.0_pr
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF ( N < 7 ) CALL ErrorExit( "OneDquadReal: too few sample points" )

      quadValue = c1*(func(1)+func(N))   + c2*(func(2)+func(N-1)) + &
                  c3*(func(3)+func(N-2)) + SUM( func(4:N-3) )

      END FUNCTION OneDquadReal

   END FUNCTION Quad2dFullgridReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad2dFullgridCmplx( TwoDfunc, x, y ) RESULT( quad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the integral of the input complex-valued function over the full
   ! two-dimensional grid, using successive 1-d quadratures. The function
   ! values are assumed to be uniformly-spaced.

   !      9 July '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: TwoDfunc
   REAL(pr),     DIMENSION(:),   INTENT(IN) :: x, y

   COMPLEX(prc) :: quad

   INTRINSIC SIZE, REAL

   ! Local variables
   COMPLEX(prc), DIMENSION(SIZE( TwoDfunc, 1 )) :: yFunc
   INTEGER  :: Nx, Ny, j
   REAL(pr) :: dx, dy
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( TwoDfunc, 2 )
   Ny = SIZE( TwoDfunc, 1 )

   dx = ( x(Nx)-x(1) )/REAL( Nx-1, KIND=pr )
   dy = ( y(Ny)-y(1) )/REAL( Ny-1, KIND=pr )

   ! First integrate along x, for each y-value...

   DO j = 1, Ny
      yFunc(j) = OneDquadCmplx( TwoDfunc(j,:), Nx )
   END DO

   ! ... then integrate along y

   quad = (dx*dy)*OneDquadCmplx( yFunc, Ny )

   CONTAINS

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FUNCTION OneDquadCmplx( func, N ) RESULT( quadValue )
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Uses a closed-interval, extended formula of order (1/N^4), from
      ! "Numerical Recipes in Fortran," Second Edition, equation (4.1.14).
      !   9 July '01
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      COMPLEX(prc), DIMENSION(:), INTENT(IN) :: func
      INTEGER,                    INTENT(IN) :: N

      COMPLEX(prc)                           :: quadValue

      INTRINSIC SUM

      ! Local variables:
      REAL(pr), PARAMETER :: c1 = 3.0_pr/8.0_pr
      REAL(pr), PARAMETER :: c2 = 7.0_pr/6.0_pr
      REAL(pr), PARAMETER :: c3 = 23.0_pr/24.0_pr
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IF ( N < 7 ) CALL ErrorExit( "OneDquadReal: too few sample points" )

      quadValue = c1*(func(1)+func(N))   + c2*(func(2)+func(N-1)) + &
                  c3*(func(3)+func(N-2)) + SUM( func(4:N-3) )

      END FUNCTION OneDquadCmplx

   END FUNCTION Quad2dFullgridCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad1DReal( vector, x, xLower, xUpper, kLower, kUpper )  &
                      RESULT( integral )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Performs 1-d quadrature using Simpson rule, with corrections for ends
   ! of the integration interval. The input vector contains function values
   ! at uniformly-spaced intervals of the independent variable, x.
   !   The matrix M results from a cubic fit to four points that straddle the
   ! integration endpoints. It is used to calculate the contributions from
   ! the partial sample intervals near xLower and xUpper.

   !      2 Nov  '00
   !     19 Jan  '01    Matrix M initialized with double-precision values
   !     16 Feb  '01    Matrix M not a PARAMETER, but re-calculated (for ABSOFT)
   !      1 Jun  '01    M (now MFIT) is private module parameter (Fortran 95)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Math_Routines

   IMPLICIT NONE

   REAL(pr), DIMENSION(:), INTENT(IN) :: vector
   REAL(pr), DIMENSION(:), INTENT(IN) :: x
   REAL(pr),               INTENT(IN) :: xLower, xUpper  ! Integration limits
   INTEGER,                INTENT(IN), OPTIONAL :: kLower, kUpper

   REAL(pr) :: integral

   INTRINSIC RESHAPE, SIZE, REAL, MODULO, SUM, DOT_PRODUCT, MATMUL

   ! Local variables
   INTEGER  :: Nx
   INTEGER  :: k, k1, k2, kStart
   REAL(pr) :: dx
   REAL(pr) :: firstPartial, firstInterval
   REAL(pr) :: lowerEnd,     upperEnd
   REAL(pr), DIMENSION(4) :: coeff1, coeff2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   dx = ( x(Nx) - x(1) )/REAL( Nx-1, KIND=pr )

   ! Find sample points k1 through k2 that lie strictly inside the
   ! integration interval.
   IF ( PRESENT( kLower ) ) THEN
      k1 = kLower
   ELSE
      k1 = FindLowerIndex( xLower, x )
   END IF
   IF ( PRESENT( kUpper ) ) THEN
      k2 = kUpper
   ELSE
      k2 = FindUpperIndex( xUpper, x )
   END IF
   IF ( k1 < 3 )    CALL ErrorExit( "Quad1DReal: not enough samples at lower end" )
   IF ( k1 > Nx-2 ) CALL ErrorExit( "Quad1DReal: not enough samples at upper end" )

   ! Integral from x(k1-1) to x(k1)
   firstInterval = dx*DOT_PRODUCT( SUM( MFIT, dim=1 ), vector(k1-2:k1+1) )

   ! Rank 1 arrays containing coefficients for curve fit to function samples
   ! at the ends of the range of integration
   coeff1 = MATMUL( MFIT, vector(k1-2:k1+1) )
   coeff2 = MATMUL( MFIT, vector(k2-1:k2+2) )

   ! Contributions from the ends of the integration interval
   ! (fractions of the sample interval)
   firstPartial = ( xLower - x(k1-1) )*EvalPoly( (xLower - x(k1-1))/dx, coeff1 )
   lowerEnd = firstInterval - firstPartial
   upperEnd = ( xUpper - x(k2) ) * EvalPoly( (xUpper - x(k2))/dx,   coeff2 )

   IF ( k2 >= k1 ) THEN

      integral = lowerEnd + upperEnd

      kStart = k1
      IF ( k2 > k1 ) THEN

         IF ( MODULO( k2-kStart, 2 ) /= 0 ) THEN
            ! Do quadrature over larger interval (to have an even number
            ! of sample intervals), then subtract the  contribution from the
            ! extra sample interval.
            kStart   = k1-1
            integral = integral - firstInterval
         END IF

         ! Use Simpson quadrature for interior points
         integral = integral +                                          &
                    dx * ( vector(kStart) + vector(k2) +                &
                        4.0_pr*SUM( vector(kStart+1:k2-1:2) ) +         &
                        2.0_pr*SUM( vector(kStart+2:k2-2:2) ) ) / 3.0_pr

      END IF

   ELSE IF ( k2 == k1-1 ) THEN
      ! Lower and upper limits lie in the same sample interval
      integral = upperEnd - firstPartial

   ELSE

      CALL ErrorExit( "Quad1DReal: indexing error" )

   END IF

   END FUNCTION Quad1DReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Quad1DCmplx( vector, x, xLower, xUpper, kLower, kUpper )  &
                       RESULT( integral )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Performs 1-d quadrature using Simpson rule, with corrections for ends
   ! of the integration interval. The input vector contains function values
   ! at uniformly-spaced intervals of the independent variable, x.
   !   The matrix M results from a cubic fit to four points that straddle the
   ! integration endpoints. It is used to calculate the contributions from
   ! the partial sample intervals near xLower and xUpper.

   !      2 Nov  '00
   !     19 Jan  '01    Matrix M initialized with double-precision values
   !     16 Feb  '01    Matrix M not a PARAMETER, but re-calculated (for ABSOFT)
   !      1 Jun  '01    M (now MFIT) is private module parameter (Fortran 95)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Math_Routines

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:), INTENT(IN) :: vector
   REAL(pr),     DIMENSION(:), INTENT(IN) :: x
   REAL(pr),                   INTENT(IN) :: xLower, xUpper ! Integration limits
   INTEGER,                    INTENT(IN), OPTIONAL :: kLower, kUpper

   COMPLEX(prc) :: integral

   INTRINSIC RESHAPE, SIZE, REAL, MODULO, SUM, DOT_PRODUCT, MATMUL

   ! Local variables
   INTEGER  :: Nx
   INTEGER  :: k, k1, k2, kStart
   REAL(pr) :: dx
   COMPLEX(prc) :: firstPartial, firstInterval
   COMPLEX(prc) :: lowerEnd,     upperEnd
   COMPLEX(prc), DIMENSION(4) :: coeff1, coeff2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   dx = ( x(Nx) - x(1) )/REAL( Nx-1, KIND=pr )

   ! Find sample points k1 through k2 that lie strictly inside the
   ! integration interval.
   IF ( PRESENT( kLower ) ) THEN
      k1 = kLower
   ELSE
      k1 = FindLowerIndex( xLower, x )
   END IF
   IF ( PRESENT( kUpper ) ) THEN
      k2 = kUpper
   ELSE
      k2 = FindUpperIndex( xUpper, x )
   END IF
   IF ( k1 < 3 )    CALL ErrorExit( "Quad1DCmplx: not enough samples at lower end" )
   IF ( k1 > Nx-2 ) CALL ErrorExit( "Quad1DCmplx: not enough samples at upper end" )

   ! Integral from x(k1-1) to x(k1)
   firstInterval = dx*DOT_PRODUCT( SUM( MFIT, dim=1 ), vector(k1-2:k1+1) )

   ! Rank 1 arrays containing coefficients for curve fit to function samples
   ! at the ends of the range of integration
   coeff1 = MATMUL( MFIT, vector(k1-2:k1+1) )
   coeff2 = MATMUL( MFIT, vector(k2-1:k2+2) )

   ! Contributions from the ends of the integration interval
   ! (fractions of the sample interval)
   firstPartial = ( xLower - x(k1-1) )*EvalPoly( (xLower - x(k1-1))/dx, coeff1 )
   lowerEnd = firstInterval - firstPartial
   upperEnd = ( xUpper - x(k2) ) * EvalPoly( (xUpper - x(k2))/dx,   coeff2 )

   IF ( k2 >= k1 ) THEN

      integral = lowerEnd + upperEnd

      kStart = k1
      IF ( k2 > k1 ) THEN

         IF ( MODULO( k2-kStart, 2 ) /= 0 ) THEN
            ! Do quadrature over larger interval (to have an even number
            ! of sample intervals), then subtract the  contribution from the
            ! extra sample interval.
            kStart   = k1-1
            integral = integral - firstInterval
         END IF

         ! Use Simpson quadrature for interior points
         integral = integral +                                          &
                    dx * ( vector(kStart) + vector(k2) +                &
                        4.0_pr*SUM( vector(kStart+1:k2-1:2) ) +         &
                        2.0_pr*SUM( vector(kStart+2:k2-2:2) ) ) / 3.0_pr

      END IF

   ELSE IF ( k2 == k1-1 ) THEN
      ! Lower and upper limits lie in the same sample interval
      integral = upperEnd - firstPartial

   ELSE

      CALL ErrorExit( "Quad1DCmplx: indexing error" )

   END IF

   END FUNCTION Quad1DCmplx



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION FindLowerIndex( xValue, xArray ) RESULT( index )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! The values in xArray should be monotonically increasing. This routine
   ! finds the index for the first xArray-value such that
   !            xArray(index) > xValue
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr)               :: xValue
   REAL(pr), DIMENSION(:) :: xArray

   INTEGER :: index

   INTRINSIC SIZE

   ! Local variables:
   INTEGER :: k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   index = 0
   DO k = 1,SIZE( xArray )
      IF ( xValue < xArray(k) )  THEN
         index = k
         EXIT
      END IF
   END DO

   IF ( index < 1 ) CALL ErrorExit( "FindLowerIndex: input data error" )

   END FUNCTION FindLowerIndex



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION FindUpperIndex( xValue, xArray ) RESULT( index )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! The values in xArray should be monotonically increasing. This routine
   ! finds the index for the largest xArray-value such that
   !            xArray(index) < xValue
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr)               :: xValue
   REAL(pr), DIMENSION(:) :: xArray

   INTEGER :: index

   INTRINSIC SIZE

   ! Local variables:
   INTEGER :: k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   index = 0
   DO k = SIZE( xArray ),1,-1
      IF ( xValue > xArray(k) ) THEN
         index = k
         EXIT
      END IF
   END DO

   IF ( index < 1 ) CALL ErrorExit( "FindUpperIndex: input data error" )

   END FUNCTION FindUpperIndex

END MODULE Quadrature_Routines
