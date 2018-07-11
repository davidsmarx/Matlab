!     Last change:  RSB   8 Oct 2002   10:48 am
MODULE Math_Routines

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds
   USE Constants
   USE Error_Exit

   IMPLICIT NONE

   PUBLIC
   PRIVATE :: EvalPolyReal,     EvalPolyCmplx
   PRIVATE :: PhaseAngleScalar, PhaseAngleVector, PhaseAngleArray

   INTERFACE EvalPoly
      MODULE PROCEDURE EvalPolyReal, &
                       EvalPolyCmplx
   END INTERFACE

   INTERFACE PhaseAngle
      MODULE PROCEDURE PhaseAngleScalar, &
      					  PhaseAngleVector, &
                       PhaseAngleArray
   END INTERFACE

 CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION EvalPolyReal( x, coeffs ) RESULT( EvalPoly )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates real-valued polynomial defined by coefficients, "coeffs"
   !
   !  value = coeffs(1) + coeffs(2)*x + coeffs(3)x**2 + ...
   !
   !  20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),               INTENT(IN) :: x
   REAL(pr), DIMENSION(:), INTENT(IN) :: coeffs
   REAL(pr)                           :: EvalPoly

   INTRINSIC SIZE

   ! Local variables
   INTEGER  :: NumCoeffs, n
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   NumCoeffs = SIZE( coeffs )
   IF ( NumCoeffs <= 0 ) THEN
      EvalPoly = 0.0_pr
    ELSE
      EvalPoly = coeffs(NumCoeffs)
      DO n = NumCoeffs-1,1,-1
         EvalPoly = x*EvalPoly + coeffs(n)
      END DO
   END IF

   END FUNCTION EvalPolyReal



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION EvalPolyCmplx( x, coeffs ) RESULT( EvalPoly )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates complex-valued polynomial defined by coefficients, "coeffs"
   !
   !  value = coeffs(1) + coeffs(2)*x + coeffs(3)x**2 + ...
   !
   !   2 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),                   INTENT(IN) :: x
   COMPLEX(prc), DIMENSION(:), INTENT(IN) :: coeffs

   COMPLEX(prc)                           :: EvalPoly

   INTRINSIC SIZE

   ! Local variables
   INTEGER  :: NumCoeffs, n
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   NumCoeffs = SIZE( coeffs )
   IF ( NumCoeffs <= 0 ) THEN
      EvalPoly = ( 0.0_pr, 0.0_pr )
    ELSE
      EvalPoly = coeffs(NumCoeffs)
      DO n = NumCoeffs-1,1,-1
         EvalPoly = x*EvalPoly + coeffs(n)
      END DO
   END IF

   END FUNCTION EvalPolyCmplx



  !=============================================================================
   FUNCTION PhaseAngleScalar( z ) RESULT( theta )
  !=============================================================================

   ! For input complex number z = abs(z)*exp(i*theta), returns theta.
   ! If REAL(z) = IMAG(z) = 0, returns theta = 0.
   !      6 June '01    R.S. Benson
   !     18 Dec  '01    Revised criterion for zero phasor

   IMPLICIT NONE

   COMPLEX(prc), INTENT(IN) :: z
   REAL(pr)                 :: theta

   INTRINSIC REAL, AIMAG, ABS, TINY, ATAN2

   ! Local variables:
   REAL(pr) :: rePart, imPart
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   rePart = REAL( z )
   imPart = AIMAG( z )

      IF ( ABS( rePart ) > TINY( 0.0_pr )  .OR.  &
      	  ABS( imPart ) > TINY( 0.0_pr )        ) THEN
         theta = ATAN2( imPart, rePart )
      ELSE
         theta = 0.0_pr     ! Zero phasor; return zero value for phase
      END IF

   END FUNCTION PhaseAngleScalar



  !=============================================================================
   FUNCTION PhaseAngleVector( z ) RESULT( theta )
  !=============================================================================

   ! For vector of complex numbers z = abs(z)*exp(i*theta),
   ! returns theta values. If REAL(z) = IMAG(z) = 0, returns theta = 0.
   !      6 June '01    R.S. Benson
   !     18 Dec  '01    Revised criterion for zero phasor

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:), INTENT(IN) :: z
   REAL(pr)    , DIMENSION(SIZE(z))       :: theta

   INTRINSIC SIZE, REAL, AIMAG, ABS, TINY, ATAN2

   ! Local variables:
   REAL(pr) :: rePart, imPart
   INTEGER  :: j
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   DO j = 1, SIZE( z )

      rePart = REAL( z(j) )
      imPart = AIMAG( z(j) )

      IF ( ABS( rePart ) > TINY( 0.0_pr )  .OR.  &
      	  ABS( imPart ) > TINY( 0.0_pr )        ) THEN
         theta(j) = ATAN2( imPart, rePart )
      ELSE
         theta(j) = 0.0_pr     ! Zero phasor; return zero value for phase
      END IF

   END DO

   END FUNCTION PhaseAngleVector



  !=============================================================================
   FUNCTION PhaseAngleArray( z ) RESULT( theta )
  !=============================================================================

   ! For 2-d array of complex numbers z = abs(z)*exp(i*theta),
   ! returns theta values. If REAL(z) = IMAG(z) = 0, returns theta = 0.
   !      6 June '01    R.S. Benson
   !     18 Dec  '01    Revised criterion for zero phasor

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN)     :: z
   REAL(pr)    , DIMENSION(SIZE(z,1),SIZE(z,2)) :: theta

   INTRINSIC SIZE, REAL, AIMAG, ABS, TINY, ATAN2

   ! Local variables:
   REAL(pr) :: rePart, imPart
   INTEGER  :: j, k
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   DO j = 1, SIZE( z, 1 )
   DO k = 1, SIZE( z, 2 )

      rePart = REAL( z(j,k) )
      imPart = AIMAG( z(j,k) )

      IF ( ABS( rePart ) > TINY( 0.0_pr )  .OR.  &
      	  ABS( imPart ) > TINY( 0.0_pr )        ) THEN
         theta(j,k) = ATAN2( imPart, rePart )
      ELSE
         theta(j,k) = 0.0_pr     ! Zero phasor; return zero value for phase
      END IF

   END DO
   END DO

   END FUNCTION PhaseAngleArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION BesselJ0( x ) RESULT( bj0 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates Bessel function J0(x).
   ! Adapted from Numerical Recipes, 2nd ed.

   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN) :: x
   REAL(pr)             :: bj0

   INTRINSIC ABS, SIN, COS, SQRT, SIGN

   ! Local variables
   REAL(pr) :: ax, xx, z
   REAL(pr) :: y
   REAL(pr), DIMENSION(5) ::  p = (/ 1.0d+0,          -0.1098628627d-2, &
                                     0.2734510407d-4, -0.2073370639d-5, &
                                     0.2093887211d-6 /)
   REAL(pr), DIMENSION(5) ::  q = (/-0.1562499995d-5,  0.1430488765d-3, &
                                    -0.6911147651d-5,  0.7621095161d-6, &
                                    -0.934945152d-7 /)
   REAL(pr), DIMENSION(6) ::  r = (/  57568490574.0d+0,   -13362590354.0d+0,  &
                                        651619640.7d+0,     -11214424.18d+0,  &
                                        77392.33017d+0,     -184.9052456d+0 /)
   REAL(pr), DIMENSION(6) ::  s = (/  57568490411.0d+0,     1029532985.0d+0,  &
                                        9494680.718d+0,      59272.64853d+0,  &
                                        267.8532712d+0,              1.0d+0 /)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( x ) < 8.0_pr ) THEN
      y = x**2
      bj0 = EvalPoly( y, r )/EvalPoly( y, s )
   ELSE
      ax = ABS( x )
      z  = 8.0_pr/ax
      y  = z**2
      xx = ax - 0.785398164_pr
      bj0 = SQRT( 0.636619772_pr/ax ) * &
                ( COS( xx )*EvalPoly( y, p ) - z*SIN( xx )*EvalPoly( y, q ) )
   END IF

   END FUNCTION BesselJ0



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION BesselJ1 ( x ) RESULT( bj1 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates Bessel function J1(x).
   ! Adapted from Numerical Recipes, 2nd ed.

   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN) :: x
   REAL(pr)             :: bj1

   INTRINSIC ABS, SIN, COS, SQRT

   ! Local variables
   REAL(pr) :: ax, xx, z
   REAL(pr) :: y
   REAL(pr), DIMENSION(5) ::  p = (/ 1.0d+0,           0.183105d-2,     &
                                    -0.3516396496d-4,  0.2457520174d-5, &
                                    -0.240337019d-6 /)
   REAL(pr), DIMENSION(5) ::  q = (/ 0.4687499995d-1,  -0.2002690873d-3, &
                                     0.8449199096d-5,  -0.88228987d-6,   &
                                     0.105787412d-6 /)
   REAL(pr), DIMENSION(6) ::  r = (/  72362614232.0d+0,    -7895059235.0d+0,  &
                                        242396853.1d+0,     -2972611.439d+0,  &
                                        15704.48260d+0,     -30.16036606d+0 /)
   REAL(pr), DIMENSION(6) ::  s = (/ 144725228442.0d+0,     2300535178.0d+0,  &
                                        18583304.74d+0,      99447.43394d+0,  &
                                        376.9991397d+0,              1.0d+0 /)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( x ) < 8.0_pr ) THEN
      y = x**2
      bj1 = x * EvalPoly( y, r )/EvalPoly( y, s )
   ELSE
      ax = ABS( x )
      z  = 8.0_pr/ax
      y  = z**2
      xx = ax - 2.35619772_pr
      bj1 = SIGN( 1.0_pr, x ) * SQRT( 0.636619772_pr / ax ) *     &
                ( COS( xx )*EvalPoly( y, p ) - z*SIN( xx )*EvalPoly( y, q ) )
   END IF

   END FUNCTION BesselJ1



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION BesselK0 ( x ) RESULT( bk0 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates Bessel function K0(x).
   ! Adapted from Numerical Recipes, 2nd ed.

   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN) :: x
   REAL(pr)             :: bk0

   INTRINSIC ABS, LOG, EXP, SQRT

   ! Local variables
   REAL(pr) :: y
   REAL(pr), DIMENSION(7) ::  p = (/-0.57721566d+0,  0.42278420d+0, &
                                     0.23069756d+0,  0.3488590d-1,  &
                                     0.262698d-2,    0.10750d-3,    &
                                     0.74d-5 /)
   REAL(pr), DIMENSION(7) ::  q = (/ 1.25331414d+0, -0.7832358d-1, &
                                     0.2189568d-1,  -0.1062446d-1, &
                                     0.587872d-2,   -0.251540d-2, &
                                     0.53208d-3 /)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( x <= 0.0_pr ) CALL ErrorExit( "BesselK0: negative input parameter" )

   IF ( x <= 2.0_pr ) THEN
      y = x*x / 4.0_pr
      bk0 = ( -LOG( x/2.0_pr ) * BesselI0( x ) ) + EvalPoly( y, p )
   ELSE
      y  = 2.0_pr / x
      bk0 =( EXP( -x ) / SQRT( x ) ) * EvalPoly( y, q )
   END IF

   END FUNCTION BesselK0



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION BesselI0 ( x ) RESULT( bi0 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates Bessel function I0(x).
   ! Adapted from Numerical Recipes, 2nd ed.

   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN) :: x
   REAL(pr)             :: bi0

   INTRINSIC ABS, EXP, SQRT

   ! Local variables
   REAL(pr) :: ax
   REAL(pr) :: y
   REAL(pr), DIMENSION(7) ::  p = (/ 1.0d+0,         3.5156229d+0,  &
                                     3.0899424d+0,   1.2067492d+0,  &
                                     0.2659732d+0,   0.360768d-1,   &
                                     0.45813d-2 /)
   REAL(pr), DIMENSION(9) ::  q = (/ 0.39894228d+0,  0.1328592d-1, &
                                     0.225319d-2,   -0.157565d-2,  &
                                     0.916281d-2,   -0.2057706d-1, &
                                     0.2635537d-1,  -0.1647633d-1, &
                                     0.392377d-2 /)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ax = ABS( x )
   IF ( ax <= 3.75_pr ) THEN
      y = x*x / 4.0_pr
      bi0 = EvalPoly( (x/3.75_pr)**2, p )
   ELSE
      bi0 =( EXP( ax ) / SQRT( ax ) ) * EvalPoly( 3.75_pr/ax, q )
   END IF

   END FUNCTION BesselI0



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Jinc( x )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates 2*J1(pi*x)/(pi*x)
   ! Adapted from Numerical Recipes, 2nd ed.

   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN) :: x
   REAL(pr)             :: Jinc

   INTRINSIC ABS, SQRT, SIN, COS

   ! Local variables
   REAL(pr) :: pix, ax, xx, z
   REAL(pr) :: y
   REAL(pr), DIMENSION(5) ::  p = (/ 1.0d+0,           0.183105d-2,     &
                                    -0.3516396496d-4,  0.2457520174d-5, &
                                    -0.240337019d-6 /)
   REAL(pr), DIMENSION(5) ::  q = (/ 0.4687499995d-1, -0.2002690873d-3, &
                                     0.8449199096d-5, -0.88228987d-6,   &
                                     0.105787412d-6 /)
   REAL(pr), DIMENSION(6) ::  r = (/  72362614232.0d+0,    -7895059235.0d+0,  &
                                        242396853.1d+0,     -2972611.439d+0,  &
                                        15704.48260d+0,     -30.16036606d+0 /)
   REAL(pr), DIMENSION(6) ::  s = (/ 144725228442.0d+0,     2300535178.0d+0,  &
                                        18583304.74d+0,      99447.43394d+0,  &
                                        376.9991397d+0,              1.0d+0 /)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   pix = PI*x
   IF ( ABS( pix ) < 8.0_pr ) THEN
      y = pix**2
      Jinc = 2.0_pr * EvalPoly( y, r )/EvalPoly( y, s )
   ELSE
      ax = ABS( pix )
      z  = 8.0_pr/ax
      y  = z**2
      xx = ax - 2.35619772_pr
      Jinc = ( 2.0_pr / pix ) * SQRT( 0.636619772_pr / ax ) *     &
             ( COS( xx )*EvalPoly( y, p ) - z*SIN( xx )*EvalPoly( y, q ) )
   END IF

   END FUNCTION Jinc



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Sinc( x )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates sin(pi*x)/(pi*x)

   !     20 Oct  '00
   !     24 May  '01    Corrected expression for "one20"
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN) :: x
   REAL(pr)             :: Sinc

   INTRINSIC SIN, ABS

   ! Local variables
   REAL(pr), PARAMETER :: SMALL_X = 0.001_pr
   REAL(pr), PARAMETER :: one6    = 1.0_pr/6.0_pr
   REAL(pr), PARAMETER :: one20   = 1.0_pr/20.0_pr
   REAL(pr)            :: pix
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   pix = PI*x

   IF ( ABS( pix ) < SMALL_X ) THEN
      Sinc = 1.0_pr - one6*pix**2 * ( 1.0_pr - one20*pix**2 )
   ELSE
      Sinc = SIN( pix ) / pix
   END IF

   END FUNCTION Sinc



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION PowerOf2( N ) RESULT( pow )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Returns the smallest exponent, pow, such that 2**pow is greater than or
   ! equal to N. If N is outside the range 1 to 2**maxPow, the routine returns
   ! a value of -1.

   !      27 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   INTEGER, INTENT(IN)  :: N
   INTEGER              :: pow

   ! Local variables
   INTEGER            :: M, twoPowM
   INTEGER, PARAMETER :: maxPow = 20
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   pow = -1
   IF ( N < 1 ) RETURN

   twoPowM =  1
   DO M = 0,maxPow
      IF ( twoPowM >= N ) THEN
         pow = M
         EXIT
      END IF
      twoPowM = 2*twoPowM
   END DO

   END FUNCTION PowerOf2



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Zernike( ZmodeValue, x, y, modeNum, ZmodeCoeff, RMS )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Evaluates a Zernike polynomial, returning an array of values.
	! Mode designations adopted from "Engineering & Laboratory Notes"
   ! from the August 1994 issue of Optics & Photonics News, by V.N. Mahajan.
   ! See also: R.J. Noll, JOSA vol.66, 207-211 (1976)
   !
   !  x             = vector of x-coordinate values, normalized to unit radius
   !  y             = vector of y-coordinate values       "      "   "    "
   !  modeNum       = mode number (Noll's ordering) between 0 and 36
   !  ZmodeCoeff    = Magnitude of this mode
   !  RMS(optional) = .TRUE. if magnitude is RMS; otherwise it is peak-to-valley
   !                  (default is peak-to-valley if RMS not present)
   !
   !  ZmodeValue(j,k) = Z[y(j),x(k)], j = 1,2,..., k = 1,2,...
   !
   !  22 May 2001
   !   9 Aug '01     RMS/peak-to-valley option; re-ordered input parameters
   !  10 Aug '01     Corrected errors in calculating rSinCos
   !  24 Aug '01     Corrected calculation for case mAz =0
   !  18 Dec '01     Subroutine, rather than function
   !   8 Oct '02     Adds calculated values to ZmodeValue array; calling routine
   !                 must do the appropriate initialization of this array
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:), INTENT(IN) :: x, y
   INTEGER,                INTENT(IN) :: modeNum
   REAL(pr),               INTENT(IN) :: ZmodeCoeff
   LOGICAL, OPTIONAL,      INTENT(IN) :: RMS

   REAL(pr), DIMENSION(SIZE(y),SIZE(x)), INTENT(IN OUT) :: ZmodeValue

   INTRINSIC ABS, EPSILON, SIZE, PRESENT, SQRT, REAL, MODULO, ALLOCATED

   ! Local variables:
   INTEGER,                     PARAMETER :: modeMax = 36
   INTEGER, DIMENSION(modeMax), PARAMETER :: n = &          ! Radial mode
      (/ 0,                                      &
      	1, 1,                                   &
      	2, 2, 2,                                &
      	3, 3, 3, 3,                             &
      	4, 4, 4, 4, 4,                          &
      	5, 5, 5, 5, 5, 5,                       &
      	6, 6, 6, 6, 6, 6, 6,                    &
      	7, 7, 7, 7, 7, 7, 7, 7 /)
   INTEGER, DIMENSION(modeMax), PARAMETER :: m = &          ! Azimuthal mode
      (/ 0,                                      &
      	1, 1,                                   &
      	0, 2, 2,                                &
      	1, 1, 3, 3,                             &
      	0, 2, 2, 4, 4,                          &
      	1, 1, 3, 3, 5, 5,                       &
      	0, 2, 2, 4, 4, 6, 6,                    &
      	1, 1, 3, 3, 5, 5, 7, 7 /)

   INTEGER  :: nRad, mAz
   INTEGER  :: j, k, mIndex
   LOGICAL  :: evenMode
   REAL(pr) :: rmCosm, rmSinm, rmCosm1, rmSinm1
   REAL(pr) :: Zcoeff, value
   REAL(pr) :: rsq
   REAL(pr), DIMENSION(:), ALLOCATABLE :: QnmCoeffs
   INTEGER                             :: nCoeffs
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( ZmodeCoeff ) <= EPSILON( 0.0_pr ) ) RETURN

   IF ( modeNum < 1  .OR.  modeNum > modeMax ) THEN
      CALL ErrorExit( "Zernike: incorrect mode specification" )
   END IF

   nRad = n(modeNum)    	! Radial mode
   mAz  = m(modeNum)    	! Azimuthal mode

   ! Normalization constant; the default assumption is that
   ! ZmodeCoeff is the peak-to-valley magnitude of the mode
   Zcoeff = ZmodeCoeff/2.0_pr

   IF ( PRESENT( RMS ) ) THEN
      IF ( RMS ) THEN
   		IF ( mAz == 0 ) THEN
      		Zcoeff = ZmodeCoeff * SQRT( REAL( nRad+1, KIND=pr ) )
   		ELSE
      		Zcoeff = ZmodeCoeff * SQRT( REAL( 2*( nRad+1 ), KIND=pr ) )
   		END IF
      END IF
   END IF

   ! Coefficients for radial polynomial
   nCoeffs = 1 + ( nRad - mAz )/2
   IF ( nCoeffs >= 2 ) THEN
   	ALLOCATE( QnmCoeffs(nCoeffs) )
   	CALL SetQnmCoeffs
   END IF

   evenMode = MODULO( modeNum, 2 ) == 0

   ! Evaluate this Zernike polynomial for all spatial points
   DO k = 1, SIZE( x )
   DO j = 1, SIZE( y )

      value = Zcoeff

      IF ( nCoeffs >= 2 ) THEN   ! Evaluate radial polynomial
	      rsq   = x(k)**2 + y(j)**2
      	value = value * EvalPoly( rsq, QnmCoeffs )
      END IF

      IF ( mAz > 0 ) THEN     	! Calculate r**m * cos/sin(m*theta)
         rmCosm = x(k)
         rmSinm = y(j)
         DO mIndex = 2,mAz
            rmCosm1 = rmCosm
            rmSinm1 = rmSinm
            rmCosm = x(k)*rmCosm1 - y(j)*rmSinm1
            rmSinm = y(j)*rmCosm1 + x(k)*rmSinm1
         END DO
         IF ( evenMode ) THEN
            value = value * rmCosm
         ELSE
            value = value * rmSinm
         END IF
      END IF

      ZmodeValue(j,k) = ZmodeValue(j,k) + value

   END DO
   END DO

   IF ( ALLOCATED( QnmCoeffs ) ) DEALLOCATE( QnmCoeffs )

   CONTAINS

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE SetQnmCoeffs
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

     !Local variables - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     !(none)
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      SELECT CASE ( nRad )
         CASE ( 2 )
            IF ( nCoeffs == 2 ) THEN    		! Defocus
               QnmCoeffs(1:nCoeffs) = (/ -1.0_pr, 2.0_pr /)
            END IF
         CASE ( 3 )
            IF ( nCoeffs == 2 ) THEN    		! Primary coma
               QnmCoeffs(1:nCoeffs) = (/ -2.0_pr, 3.0_pr /)
            END IF
         CASE ( 4 )
            IF ( nCoeffs == 2 ) THEN    		! Secondary astigmatism
               QnmCoeffs(1:nCoeffs) = (/ -3.0_pr, 4.0_pr /)
            ELSE IF ( nCoeffs == 3 ) THEN    ! Primary spherical
               QnmCoeffs(1:nCoeffs) = (/ 1.0_pr, -6.0_pr, 6.0_pr /)
            END IF
         CASE ( 5 )
            IF ( nCoeffs == 2 ) THEN
               QnmCoeffs(1:nCoeffs) = (/ -4.0_pr, 5.0_pr /)
            ELSE IF ( nCoeffs == 3 ) THEN    ! Secondary coma
               QnmCoeffs(1:nCoeffs) = (/ 3.0_pr, -12.0_pr, 10.0_pr /)
            END IF
         CASE ( 6 )
            IF ( nCoeffs == 2 ) THEN
               QnmCoeffs(1:nCoeffs) = (/ -5.0_pr, 6.0_pr /)
            ELSE IF ( nCoeffs == 3 ) THEN
               QnmCoeffs(1:nCoeffs) = (/ 6.0_pr, -20.0_pr, 15.0_pr /)
            ELSE IF ( nCoeffs == 4 ) THEN    ! Secondary spherical
               QnmCoeffs(1:nCoeffs) = (/ -1.0_pr, 12.0_pr, -30.0_pr, 20.0_pr /)
            END IF
         CASE ( 7 )
            IF ( nCoeffs == 2 ) THEN
               QnmCoeffs(1:nCoeffs) = (/ -6.0_pr, 7.0_pr /)
            ELSE IF ( nCoeffs == 3 ) THEN
               QnmCoeffs(1:nCoeffs) = (/ 10.0_pr, -30.0_pr, 21.0_pr /)
            ELSE IF ( nCoeffs == 4 ) THEN    ! Tertiary coma
               QnmCoeffs(1:nCoeffs) = (/ -4.0_pr, 30.0_pr, -60.0_pr, 35.0_pr /)
            END IF
         CASE DEFAULT
            CALL ErrorExit( "SetQnmCoeffs: incorrect radial mode specification" )
      END SELECT

      END SUBROUTINE SetQnmCoeffs

	END SUBROUTINE Zernike

END MODULE Math_Routines
