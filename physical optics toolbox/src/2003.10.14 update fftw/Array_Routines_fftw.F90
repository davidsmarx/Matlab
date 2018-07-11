!     Last change:  RSB  23 Sep 2002    4:01 pm
MODULE Array_Routines

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds
   USE Constants

   IMPLICIT NONE

   PUBLIC

   PRIVATE :: MakeRealArray,              MakeComplexArray
   PRIVATE :: PadRealArray,               PadComplexArray
   PRIVATE :: UnPadRealArray,             UnPadComplexArray

   PRIVATE :: XYZ_RotationMatrix
   PRIVATE :: ReflectRealAboutXaxis,      ReflectRealAboutYaxis
   PRIVATE :: ReflectComplexAboutXaxis,   ReflectComplexAboutYaxis

   INTERFACE MakeArray
      MODULE PROCEDURE MakeRealArray, MakeComplexArray
   END INTERFACE

   INTERFACE PadArray
      MODULE PROCEDURE PadRealArray, PadComplexArray
   END INTERFACE

   INTERFACE UnPadArray
      MODULE PROCEDURE UnPadRealArray, UnPadComplexArray
   END INTERFACE

   INTERFACE RotationMatrix
      MODULE PROCEDURE XYZ_RotationMatrix, General_RotationMatrix
   END INTERFACE

   INTERFACE ReflectAboutXaxis
      MODULE PROCEDURE ReflectRealAboutXaxis, ReflectComplexAboutXaxis
   END INTERFACE

   INTERFACE ReflectAboutYaxis
      MODULE PROCEDURE ReflectRealAboutYaxis, ReflectComplexAboutYaxis
   END INTERFACE

   INTEGER, DIMENSION(2), PARAMETER :: rowOrder       = (/ 2, 1 /)
   INTEGER, DIMENSION(2), PARAMETER :: two_by_two     = (/ 2, 2 /)
   INTEGER, DIMENSION(2), PARAMETER :: three_by_three = (/ 3, 3 /)

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION CoordinateVector( N, dx, xc ) RESULT( vector )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! CoordinateVector: returns 1-d array of evenly-spaced coordinate values

   !     N  = Number of values desired
   !     dx = interval between values
   !     xc = (OPTIONAL) coordinate value at center grid point

   !     17 Dec '01     Added optional coordinate value at "center" (i.e., Nc)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   INTEGER,  INTENT(IN)             :: N
   REAL(pr), INTENT(IN)             :: dx
   REAL(pr), INTENT(IN), OPTIONAL   :: xc
   REAL(pr), DIMENSION(N)           :: vector

   INTRINSIC REAL, PRESENT

   !  Local variables
   INTEGER k, Nc
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nc = ( N/2 ) + 1        ! Coordinate center

   DO k = 1,N
      vector(k) = dx * REAL( (k-Nc), KIND=pr )
   END DO
   IF ( PRESENT( xc ) ) vector = vector + xc

   END FUNCTION CoordinateVector



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION MakeRealArray( rows, cols, realVector ) RESULT( realArray )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! MakeRealArray: reshapes input vector into a 2-d array, i.e.,
   !     realArray(i,j) = realVector([i-1]*cols+j), i = 1,rows, j = 1,cols

   !     25 May  '01
   ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   USE Error_Exit

   IMPLICIT NONE

   INTEGER,                INTENT(IN) :: rows, cols
   REAL(pr), DIMENSION(:), INTENT(IN) :: realVector

   REAL(pr), DIMENSION(rows,cols) :: realArray

	! Local variables
   ! (none)
   ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   IF ( rows*cols == SIZE( realVector ) ) THEN
      realArray = RESHAPE( realVector, (/ rows, cols /), ORDER=rowOrder )
   ELSE
      CALL ErrorExit( "MakeRealArray: required array size doesn't match input" )
   END IF

   END FUNCTION MakeRealArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION MakeComplexArray( rows, cols, cmplxVector ) RESULT( cmplxArray )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! MakeComplexArray: reshapes input vector into a 2-d array, i.e.,
   !     cmplxArray(i,j) = cmplxVector([i-1]*cols+j), i = 1,rows, j = 1,cols

   !     25 May  '01
   ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   USE Error_Exit

   IMPLICIT NONE

   INTEGER,                    INTENT(IN) :: rows, cols
   COMPLEX(prc), DIMENSION(:), INTENT(IN) :: cmplxVector

   COMPLEX(prc), DIMENSION(rows,cols) :: cmplxArray

	! Local variables
   ! (none)
   ! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

   IF ( rows*cols == SIZE( cmplxVector ) ) THEN
      cmplxArray = RESHAPE( cmplxVector, (/ rows, cols /), ORDER=rowOrder )
   ELSE
      CALL ErrorExit( "MakeComplexArray: required array size doesn't match input" )
   END IF

   END FUNCTION MakeComplexArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	FUNCTION PadRealArray( array, Nrows, Ncolumns, rowOffset, colOffset ) &
				RESULT( paddedArray )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	! PadRealArray: adds zero-value samples to an input array. The center
	! of the input array is placed at the center of the new array if no
	! offsets are specified.

	!		24 July '00		Added optional offset inputs
	!		12 Mar  '01		Fortran 90 version

	!    INPUT:
	!           array     = original 2-d, real-valued array
	!           Nrows     = final number of rows of returned array
	!           Ncolumns  = final number of columns
	!optional	rowOffset = row offset from center of input array that will
	!								become the center of the unpadded array
	!optional	colOffset = column offset from center of input array that will
	!								become the center of the unpadded array
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	USE Error_Exit

	IMPLICIT NONE

	REAL(pr), DIMENSION(:,:), INTENT(IN) :: array
	INTEGER,                  INTENT(IN) :: Nrows,     Ncolumns
	INTEGER,  OPTIONAL,       INTENT(IN) :: rowOffset, colOffset

	REAL(pr), DIMENSION(Nrows,Ncolumns) :: paddedArray

	INTRINSIC SIZE, PRESENT, MAX, MIN

	! Local variables
	INTEGER :: M, N, mdiff, ndiff, mstart, mstop, nstart, nstop
   INTEGER :: mIndex, nIndex
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		M = SIZE( array, 1 )
		N = SIZE( array, 2 )

		IF ( Nrows < M  .OR.  Ncolumns < N ) THEN
			CALL ErrorExit( "Improper input for PadRealArray" )
		ELSE IF ( Nrows == M  .AND.  Ncolumns == N ) THEN
			paddedArray = array
		ELSE
			mdiff = (Nrows/2)    - (M/2)
			ndiff = (Ncolumns/2) - (N/2)

			IF ( PRESENT( rowOffset ) .AND. PRESENT( colOffset ) ) THEN
				mdiff = mdiff + rowOffset
				ndiff = ndiff + colOffset
			END IF

			! Place old array in the center of the new array. In the case of
			! non-zero offsets, points of the input array are ignored if they
			! fall outside the new array.

	      paddedArray(1:Nrows,1:Ncolumns) = 0.0_pr
			mstart = MAX( 1 , 1-mdiff )
			mstop  = MIN( M , Nrows-mdiff )
			nstart = MAX( 1 , 1-ndiff )
			nstop  = MIN( N , Ncolumns-ndiff )

			DO mIndex = mstart,mstop
				DO nIndex = nstart,nstop
					paddedArray(mIndex+mdiff,nIndex+ndiff) = array(mIndex,nIndex)
				END DO
			END DO

		END IF

	END FUNCTION PadRealArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	FUNCTION PadComplexArray( array, Nrows, Ncolumns, rowOffset, colOffset ) &
				RESULT( paddedArray )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	! PadRealArray: adds zero-value samples to an input array. The center
	! of the input array is placed at the center of the new array if no
	! offsets are specified.

	!		24 July '00		Added optional offset inputs
	!		12 Mar  '01		Fortran 90 version

	!    INPUT:
	!           array     = original 2-d, real-valued array
	!           Nrows     = final number of rows of returned array
	!           Ncolumns  = final number of columns
	!optional	rowOffset = row offset from center of input array that will
	!								become the center of the unpadded array
	!optional	colOffset = column offset from center of input array that will
	!								become the center of the unpadded array
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	USE Error_Exit

	IMPLICIT NONE

	COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: array
	INTEGER,                      INTENT(IN) :: Nrows,     Ncolumns
	INTEGER,      OPTIONAL,       INTENT(IN) :: rowOffset, colOffset

	COMPLEX(prc), DIMENSION(Nrows,Ncolumns) :: paddedArray

	INTRINSIC SIZE, PRESENT, MAX, MIN

	! Local variables
	INTEGER :: M, N, mdiff, ndiff, mstart, mstop, nstart, nstop
   INTEGER :: mIndex, nIndex
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		M = SIZE( array, 1 )
		N = SIZE( array, 2 )

		IF ( Nrows < M  .OR.  Ncolumns < N ) THEN
			CALL ErrorExit( "Improper input for PadComplexArray" )
		ELSE IF ( Nrows == M  .AND.  Ncolumns == N ) THEN
			paddedArray = array
		ELSE
			mdiff = (Nrows/2)    - (M/2)
			ndiff = (Ncolumns/2) - (N/2)

			IF ( PRESENT( rowOffset ) .AND. PRESENT( colOffset ) ) THEN
				mdiff = mdiff + rowOffset
				ndiff = ndiff + colOffset
			END IF

			! Place old array in the center of the new array. In the case of
			! non-zero offsets, points of the input array are ignored if they
			! fall outside the new array.

	      paddedArray(1:Nrows,1:Ncolumns) = ( 0.0_pr, 0.0_pr )
			mstart = MAX( 1 , 1-mdiff )
			mstop  = MIN( M , Nrows-mdiff )
			nstart = MAX( 1 , 1-ndiff )
			nstop  = MIN( N , Ncolumns-ndiff )

			DO mIndex = mstart,mstop
				DO nIndex = nstart,nstop
					paddedArray(mIndex+mdiff,nIndex+ndiff) = array(mIndex,nIndex)
				END DO
			END DO

		END IF

	END FUNCTION PadComplexArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	FUNCTION UnPadRealArray( array, Nrows, Ncolumns, rowOffset, colOffset ) &
				RESULT( unPaddedArray )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	! UnPadRealArray: Eliminates "padding" around edges of the input array

	!		24 July '00		Added optional offset inputs
	!		12 Mar  '01		Fortran 90 version

	!    INPUT:
	!           array     = original 2-d, real-valued array
	!           Nrows     = final number of rows of returned array
	!           Ncolumns  = final number of columns
	!optional	rowOffset = row offset from center of input array that will
	!								become the center of the unpadded array
	!optional	colOffset = column offset from center of input array that will
	!								become the center of the unpadded array
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	USE Error_Exit

	IMPLICIT NONE

	REAL(pr), DIMENSION(:,:), INTENT(IN) :: array
	INTEGER,                  INTENT(IN) :: Nrows,     Ncolumns
	INTEGER,  OPTIONAL,       INTENT(IN) :: rowOffset, colOffset

	REAL(pr), DIMENSION(Nrows,Ncolumns) :: unPaddedArray

	INTRINSIC SIZE, PRESENT, MAX, MIN

	! Local variables
	INTEGER :: M, N, mdiff, ndiff, mstart, mstop, nstart, nstop
   INTEGER :: mIndex, nIndex
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		M = SIZE( array, 1 )
		N = SIZE( array, 2 )

		IF ( Nrows > M  .OR.  Ncolumns > N ) THEN
			CALL ErrorExit( "Improper input for UnPadRealArray" )
		ELSE IF ( Nrows == M  .AND.  Ncolumns == N ) THEN
			unPaddedArray = array
		ELSE
			mdiff = (M/2) - (Nrows/2)
			ndiff = (N/2) - (Ncolumns/2)

			IF ( PRESENT( rowOffset ) .AND. PRESENT( colOffset ) ) THEN
				mdiff = mdiff + rowOffset
				ndiff = ndiff + colOffset
			END IF

			! Select the new array from the "center" of the old array, using zero
			! values for points that fall outside the input array

	      unPaddedArray(1:Nrows,1:Ncolumns) = 0.0_pr
			mstart = MAX( 1 , 1-mdiff )
			mstop  = MIN( Nrows , M-mdiff )
			nstart = MAX( 1 , 1-ndiff )
			nstop  = MIN( Ncolumns , N-ndiff )

			DO mIndex = mstart,mstop
				DO nIndex = nstart,nstop
					unPaddedArray(mIndex,nIndex) = array(mIndex+mdiff,nIndex+ndiff)
				END DO
			END DO

		END IF

	END FUNCTION UnPadRealArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
	FUNCTION UnPadComplexArray( array, Nrows, Ncolumns, rowOffset, colOffset ) &
				RESULT( unPaddedArray )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

	! UnPadComplexArray: Eliminates "padding" around edges of the input array

	!		24 July '00		Added optional offset inputs
	!		12 Mar  '01		Fortran 90 version

	!    INPUT:
	!           array     = original 2-d, real-valued array
	!           Nrows     = final number of rows of returned array
	!           Ncolumns  = final number of columns
	!optional	rowOffset = row offset from center of input array that will
	!								become the center of the unpadded array
	!optional	colOffset = column offset from center of input array that will
	!								become the center of the unpadded array
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

	USE Error_Exit

	IMPLICIT NONE

	COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: array
	INTEGER,                      INTENT(IN) :: Nrows,     Ncolumns
	INTEGER,      OPTIONAL,       INTENT(IN) :: rowOffset, colOffset

	COMPLEX(prc), DIMENSION(Nrows,Ncolumns) :: unPaddedArray

	INTRINSIC SIZE, PRESENT, MAX, MIN

	! Local variables
	INTEGER :: M, N, mdiff, ndiff, mstart, mstop, nstart, nstop
   INTEGER :: mIndex, nIndex
	! . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .

		M = SIZE( array, 1 )
		N = SIZE( array, 2 )

		IF ( Nrows > M  .OR.  Ncolumns > N ) THEN
			CALL ErrorExit( "Improper input for UnPadComplexArray" )
		ELSE IF ( Nrows == M  .AND.  Ncolumns == N ) THEN
			unPaddedArray = array
		ELSE
			mdiff = (M/2) - (Nrows/2)
			ndiff = (N/2) - (Ncolumns/2)

			IF ( PRESENT( rowOffset ) .AND. PRESENT( colOffset ) ) THEN
				mdiff = mdiff + rowOffset
				ndiff = ndiff + colOffset
			END IF

			! Select the new array from the "center" of the old array, using zero
			! values for points that fall outside the input array

	      unPaddedArray(1:Nrows,1:Ncolumns) = ( 0.0_pr, 0.0_pr )
			mstart = MAX( 1 , 1-mdiff )
			mstop  = MIN( Nrows , M-mdiff )
			nstart = MAX( 1 , 1-ndiff )
			nstop  = MIN( Ncolumns , N-ndiff )

			DO mIndex = mstart,mstop
				DO nIndex = nstart,nstop
					unPaddedArray(mIndex,nIndex) = array(mIndex+mdiff,nIndex+ndiff)
				END DO
			END DO

		END IF

	END FUNCTION UnPadComplexArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE QuadraticPhase( A, x, y, const, cxLin, cyLin, cQuad )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Applies constant, linear and quadratic phase to input complex-valued array:
   !
   !  phase(x,y) = const + x*cxLin + y*cyLin + cQuad*(x**2 + y**2)
   !
   !  A(output) = A(input)*exp(i*phase)
   !
   !  x is a vector of sample coordinates for the columns of A;
   !  y is a vector "    "        "        "   "  rows of A.
   !
   !     16 July '01
   !     30 July '01    Corrected expression for yPhasor(:); added j-loop index
   !      6 Aug  '01    Corrected testing of piston term ("const")
   !     13 Dec  '01    Revised testing for negligible phase (bug fix)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Error_Exit

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: A
   REAL(pr),     DIMENSION(:),   INTENT(IN)     :: x, y
   REAL(pr),                     INTENT(IN)     :: const, cxLin, cyLin, cQuad

   INTRINSIC SIZE, ABS, EPSILON, EXP

   ! Local variables:
   COMPLEX(prc)                       :: xPhasor
   COMPLEX(prc), DIMENSION(SIZE( y )) :: yPhasor
   INTEGER                            :: Ny, Nx, j, k
   REAL(pr)                           :: xMax, yMax, crit
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx = SIZE( x )
   Ny = SIZE( y )

   IF ( SIZE( A, 2 ) /= Nx  .OR.  SIZE( A, 1 ) /= Ny ) THEN
      CALL ErrorExit( "QuadraticPhase: dimensions do not match" )
   END IF

   xMax = ABS( x(Nx) - x(1) )
   yMax = ABS( y(Ny) - y(1) )
   crit = EPSILON(xMax)
   IF ( ABS( const ) < crit                      .AND.  &
   	  ABS(cxLin*xMax) + ABS(cyLin*yMax) < crit .AND.  &
		  ABS(cQuad*(xMax**2 + yMax**2))    < crit        ) RETURN

   yPhasor = EXP( i*( const + y*( cyLin + cQuad*y ) ) )

   DO k = 1,Nx
   	xPhasor = EXP( i * x(k)*( cxLin + cQuad*x(k) ) )
   	DO j = 1,Ny
   		A(j,k)  = A(j,k) * ( xPhasor*yPhasor(j) )
   	END DO
   END DO

   END SUBROUTINE QuadraticPhase



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION R_2d( angle ) RESULT( rotationMatrix )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Two-dimensional rotation matrix
   !  rotates an object counter-clockwise
   !  by the input angle (radians).
   !     30 July '98
   !     11 Oct  '00    Fortran 95 version
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), INTENT(IN)     :: angle
   REAL(pr), DIMENSION(2,2) :: rotationMatrix

   INTRINSIC RESHAPE

   !  Local variables
   REAL(pr) cosine, sine
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cosine = cos( angle )
   sine   = sin( angle )

   rotationMatrix = RESHAPE( (/ cosine,  -sine, &
                                 sine,   cosine /), two_by_two, &
                                                    ORDER = rowOrder )
   END FUNCTION R_2d



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION XYZ_RotationMatrix( angle, axis ) RESULT( matrix )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Retruns a rotation matrix for rotating a 3-d vector counter-clockwise
   ! about the specified axis: "x-axis" "y-axis" or "z-axis"
   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),         INTENT(IN) :: angle     ! Rotation angle (radians)
   CHARACTER(LEN=6), INTENT(IN) :: axis
   REAL(pr), DIMENSION(3,3)     :: matrix

   INTRINSIC SIN, COS, RESHAPE

   ! Local variables
   REAL(pr) :: cosine, sine
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cosine = COS( angle )
   sine   = SIN( angle )

   SELECT CASE ( axis )
      CASE ( "x-axis" )
         matrix = RESHAPE( (/ 1.0_dp,  0.0_dp,  0.0_dp,  &
                              0.0_dp,  cosine,  -sine,   &
                              0.0_dp,   sine,   cosine /), three_by_three, &
                                                           ORDER = rowOrder )
      CASE ( "y-axis" )
         matrix = RESHAPE( (/ cosine,  0.0_dp,   sine,  &
                              0.0_dp,  1.0_dp,  0.0_dp,   &
                              -sine,   0.0_dp,  cosine /), three_by_three, &
                                                           ORDER = rowOrder )
      CASE ( "z-axis" )
         matrix = RESHAPE( (/ cosine,  -sine,  0.0_dp,  &
                               sine,   cosine, 0.0_dp,   &
                              0.0_dp,  0.0_dp, 1.0_dp /), three_by_three, &
                                                           ORDER = rowOrder )
      CASE DEFAULT
   END SELECT

   END FUNCTION XYZ_RotationMatrix



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION General_RotationMatrix( angle, axis ) RESULT( matrix )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Retruns an orthogonal matrix for rotating a 3-d vector counter-clockwise
   ! about an axis specified by an input vector. For reference, see
   !   E.J. Konopinski, "Classical Descriptions of Motion," W.H. Freeman (1969),
   !   sec. 9.5, and I.N. Bronshtein & K.A. Semendyayev, "Handbook of
   !   Mathematics," 3rd ed., Van Nostrand Reinhold, N.Y. (1985), p. 196
   !     19 Sept '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),               INTENT(IN) :: angle     ! Rotation angle (radians)
   REAL(pr), DIMENSION(3), INTENT(IN) :: axis      ! Vector defining axis
   REAL(pr), DIMENSION(3,3)           :: matrix

   INTRINSIC SIN, COS, SQRT, DOT_PRODUCT

   ! Local variables
   REAL(pr) :: cosine, sine
   REAL(pr) :: absAxis
   REAL(pr) :: alpha, beta, gamma
   REAL(pr) :: temp
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cosine = COS( angle )
   sine   = SIN( angle )

   absAxis = SQRT( DOT_PRODUCT( axis, axis ) )
   alpha = axis(1)/absAxis
   beta  = axis(2)/absAxis
   gamma = axis(3)/absAxis

   temp = alpha*(1.0_pr - cosine)
   matrix(1,1) =   cosine   + alpha*temp
   matrix(1,2) = gamma*sine + beta*temp
   matrix(1,3) = -beta*sine + gamma*temp

   temp = beta*(1.0_pr - cosine)
   matrix(2,1) = -gamma*sine + alpha*temp
   matrix(2,2) =    cosine   + beta*temp
   matrix(2,3) =  alpha*sine + gamma*temp

   temp = gamma*(1.0_pr - cosine)
   matrix(3,1) =   beta*sine + alpha*temp
   matrix(3,2) = -alpha*sine + beta*temp
   matrix(3,3) =    cosine   + gamma*temp

   END FUNCTION General_RotationMatrix



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION X_Rotation( angle ) RESULT( matrix )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Retruns a rotation matrix for rotating a 3-d vector counter-clockwise
   ! about the x-axis.
   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),             INTENT(IN) :: angle    ! Rotation angle (radians)
   REAL(pr), DIMENSION(3,3)         :: matrix

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   matrix = XYZ_RotationMatrix( angle, "x-axis" )

   END FUNCTION X_Rotation



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Y_Rotation( angle ) RESULT( matrix )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Retruns a rotation matrix for rotating a 3-d vector counter-clockwise
   ! about the y-axis.
   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),             INTENT(IN) :: angle    ! Rotation angle (radians)
   REAL(pr), DIMENSION(3,3)         :: matrix

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   matrix = XYZ_RotationMatrix( angle, "y-axis" )

   END FUNCTION Y_Rotation



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION Z_Rotation( angle ) RESULT( matrix )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Retruns a rotation matrix for rotating a 3-d vector counter-clockwise
   ! about the z-axis.
   !     20 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),             INTENT(IN) :: angle    ! Rotation angle (radians)
   REAL(pr), DIMENSION(3,3)         :: matrix

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   matrix = XYZ_RotationMatrix( angle, "z-axis" )

   END FUNCTION Z_Rotation



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ReflectRealAboutXaxis( array )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Input array represents real-valued f(y,x); ouput is f(-y,x). The coordinate
   !  origin is the grid center at (ny/2)+1, (nx/2)+1.
   !      2 Feb '01     Robert Benson
   !      3 Dec '01     Explicit ncol-loop
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array

   INTRINSIC SIZE

   ! Local variables
   INTEGER  :: Ny, Ncenter, nrow
   INTEGER  :: Nx, ncol
   REAL(pr) :: temp
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx      = SIZE( array, 2 )
   Ny      = SIZE( array, 1 )
   Ncenter = 1 + (Ny/2)

   ! Swap rows on either side of the "center" row
   DO nrow = 1, Ny-Ncenter
   DO ncol = 1, Nx
      temp = array(Ncenter+nrow,ncol)
      array(Ncenter+nrow,ncol) = array(Ncenter-nrow,ncol)
      array(Ncenter-nrow,ncol) = temp
   END DO
   END DO

   END SUBROUTINE ReflectRealAboutXaxis


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ReflectComplexAboutXaxis( array )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Input array represents complex-valued f(y,x); ouput is f(-y,x).
   !  The coordinate origin is the grid center at (ny/2)+1, (nx/2)+1.
   !      2 Feb '01     Robert Benson
   !      3 Dec '01     Explicit ncol-loop
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array

   INTRINSIC SIZE

   ! Local variables
   INTEGER      :: Ny, Ncenter, nrow
   INTEGER      :: Nx, ncol
   COMPLEX(prc) :: temp
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nx      = SIZE( array, 2 )
   Ny      = SIZE( array, 1 )
   Ncenter = 1 + (Ny/2)

   ! Swap rows on either side of the "center" row
   DO nrow = 1, Ny-Ncenter
   DO ncol = 1, Nx
      temp = array(Ncenter+nrow,ncol)
      array(Ncenter+nrow,ncol) = array(Ncenter-nrow,ncol)
      array(Ncenter-nrow,ncol) = temp
   END DO
   END DO

   END SUBROUTINE ReflectComplexAboutXaxis




  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ReflectRealAboutYaxis( array )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Input array represents real-valued f(y,x); ouput is f(y,-x). The coordinate
   !  origin is the grid center at (ny/2)+1, (nx/2)+1.
   !      2 Feb '01     Robert Benson
   !      3 Dec '01     Explicit nrow-loop
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr), DIMENSION(:,:), INTENT(IN OUT) :: array

   INTRINSIC SIZE

   ! Local variables
   INTEGER  :: Nx, Ncenter, ncol
   INTEGER  :: Ny, nrow
   real(pr) :: temp
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny      = SIZE( array, 1 )
   Nx      = SIZE( array, 2 )
   Ncenter = 1 + (Nx/2)

   ! Swap columns on either side of the "center" column
   DO ncol = 1, Nx-Ncenter
   DO nrow = 1, Ny
      temp = array(nrow,Ncenter+ncol)
      array(nrow,Ncenter+ncol) = array(nrow,Ncenter-ncol)
      array(nrow,Ncenter-ncol) = temp
   END DO
   END DO

   END SUBROUTINE ReflectRealAboutYaxis



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ReflectComplexAboutYaxis( array )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Input array represents complex-valued f(y,x); ouput is f(y,-x).
   !  The coordinate origin is the grid center at (ny/2)+1, (nx/2)+1.
   !      2 Feb '01     Robert Benson
   !      3 Dec '01     Explicit nrow-loop
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array

   INTRINSIC SIZE

   ! Local variables
   INTEGER      :: Nx, Ncenter, ncol
   INTEGER      :: Ny, nrow
   COMPLEX(prc) :: temp
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny      = SIZE( array, 1 )
   Nx      = SIZE( array, 2 )
   Ncenter = 1 + (Nx/2)

   ! Swap columns on either side of the "center" column
   DO ncol = 1, Nx-Ncenter
   DO nrow = 1, Ny
      temp = array(nrow,Ncenter+ncol)
      array(nrow,Ncenter+ncol) = array(nrow,Ncenter-ncol)
      array(nrow,Ncenter-ncol) = temp
   END DO
   END DO

   END SUBROUTINE ReflectComplexAboutYaxis



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ReflectThroughOrigin( array, xOrigin, yOrigin )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Input array represents f(y,x); ouput is f(-y,-x). The coordinate origin is
   !  relative to the grid center at (ny/2)+1, (nx/2)+1. The algorithm uses
   !  theorems from Fourier analysis:
   !     If F{} denotes the Fourier transform, then F{F{f(y,x)}} = f(-y,-x),
   !     and F{f(y-yc,x-xc)} = exp(-2*pi*i*[xFreq*xc+yFreq*yc])*F{f(y,x)}

   !     27 Apr '00     Robert Benson
   !      3 Nov '00     Fortran 90/95
   !     14 Nov '00     Assumes FFT (or its inverse) is NOT normalized by 1/N
   !      2 Feb '01     xOrigin and yOrigin in units of sample spacing;
   !                    calls ReflectAboutXaxis and ReflectAboutYaxis
   !     30 Nov '01     Use DO loops for multiplying array by constant
   !     18 Dec '01     Calls revised FFTshift

   ! INPUT
   !  array             2-d array, f(y,x)
   !  xOrigin, yOrigin  Location of reflection point in units of sample spacing
   !                    -- OPTIONAL INPUT --
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Discrete_Transforms
   USE Error_Exit

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT)       :: array
   REAL(pr),                     INTENT(IN), OPTIONAL :: xOrigin, yOrigin

   INTRINSIC SIZE, PRESENT, EPSILON, EXP, CONJG

   ! Local variables
   INTEGER  :: ny, nx, j, k
   REAL(pr) :: const
   REAL(pr) :: xo, yo
   REAL(pr) :: xo_df, yo_df
   COMPLEX(prc), DIMENSION(SIZE(array,1)) :: yPhasor
   COMPLEX(prc), DIMENSION(SIZE(array,2)) :: xPhasor
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ny = SIZE( array, 1 )
   nx = SIZE( array, 2 )
   const = 1.0_pr/REAL( nx*ny, KIND=pr)

   xo = 0.0_pr
   yo = 0.0_pr
   IF ( PRESENT( xOrigin ) ) xo = xOrigin
   IF ( PRESENT( yOrigin ) ) yo = yOrigin

   IF ( abs( xo ) > EPSILON(xo)  .OR.  abs( yo ) > EPSILON(yo)  ) THEN

      CALL FFTshift( array, inverse=.TRUE. )
      CALL FFT2D( array )

      ! Shift f(y,x) so origin coincides with grid center, reflect, then shift back
      ! Shifting is done in frequency space by multiplying the Fourier transform
      ! of f(y,x) by an array of complex, unit-magnitude phasors.

      IF ( abs( xo ) > EPSILON(xo) ) THEN
         xo_df = xo/REAL( nx, KIND=pr)
         xPhasor = EXP( i*( TWOPI*CoordinateVector( nx, xo_df ) ) )
         CALL FFTshift( xPhasor, inverse=.TRUE. )
         DO k = 1,nx
         	DO j = 1,ny
               array(j,k) = xPhasor(k)*array(j,k)
            END DO
         END DO
      END IF
      IF ( abs( yo ) > EPSILON(yo) ) THEN
         yo_df = yo/REAL( nx, KIND=pr)
         yPhasor = EXP( i*( TWOPI*CoordinateVector( ny, yo_df ) ) )
         CALL FFTshift( yPhasor, inverse=.TRUE. )
         DO k = 1,nx
         	DO j = 1,ny
               array(j,k) = yPhasor(j)*array(j,k)
            END DO
         END DO
      END IF
      CALL FFT2D( array )

      ! Shift back to initial origin

      CALL FFT2D( array )
      IF ( abs( xo ) > EPSILON(xo) ) THEN
         xPhasor = CONJG( xPhasor )
         DO k = 1,nx
         	DO j = 1,ny
               array(j,k) = xPhasor(k)*array(j,k)
            END DO
         END DO
      END IF
      IF ( abs( yo ) > EPSILON(yo) ) THEN
         yPhasor = CONJG( yPhasor )
         DO k = 1,nx
         	DO j = 1,ny
               array(j,k) = yPhasor(j)*array(j,k)
            END DO
         END DO
      END IF
      CALL FFT2D( array, inverse=.TRUE. )

      DO k = 1,nx
      	DO j = 1,ny
      		array(j,k) = (const**2)*array(j,k)
         END DO
      END DO
      CALL FFTshift( array, inverse=.FALSE. )

   ELSE

      CALL ReflectComplexAboutXaxis( array )
      CALL ReflectComplexAboutYaxis( array )

   END IF

   END SUBROUTINE ReflectThroughOrigin


!DEC$ IF (.FALSE.)

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE TrigInterp2D( g,    dx,     dy,     x1,     y1,    &
                            gInt, deltaX, deltaY, xFirst, yFirst )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Two-dimensional trigonometric interpolation. A chirp-z algorithm is applied
   !  to the Fourier spectrum of the input function, which is shifted (using the
   !  Shift Theorem) to properly center it in the interpolation interval.
   !     21 Jun '00		R.S. Benson
   !      3 Nov '00		Fortran 90/95
   !     30 Nov '01     Use xConst and yConst to normalize FFT2D
   !     18 Dec '01     Calls revised FFTshift

   !  g                 g(y,x): array of input values (uniform spacing in x and y)
   !  dx, dy            Sample spacing for g
   !  x1, y1            Coordinate values for first sample of g
   !  gInt              Output array of interpolated values
   !  deltaX, deltaY    Sample spacing for interpolated output
   !  xFirst, yFirst    Coordinate values at beginning of interpolation interval
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Discrete_Transforms

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN) :: g
   REAL(pr),                     INTENT(IN) :: dx, dy
   REAL(pr),                     INTENT(IN) :: x1, y1

   COMPLEX(prc), DIMENSION(:,:), INTENT(OUT) :: gInt
   REAL(pr),                     INTENT(IN)  :: deltaX, deltaY
   REAL(pr),                     INTENT(IN)  :: xFirst, yFirst

   INTRINSIC SIZE, REAL, EXP

   ! Local variables
   INTEGER  :: nGy, nGx, Nx, Ny
   REAL(pr) :: xCenter, yCenter
   REAL(pr) :: xShift,  yShift
   REAL(pr) :: xConst,  yConst
   REAL(pr) :: dfx,     dfy
   INTEGER  :: j, k
   COMPLEX(prc), DIMENSION(:),   ALLOCATABLE :: yPhasor
   COMPLEX(prc), DIMENSION(:),   ALLOCATABLE :: xPhasor
   COMPLEX(prc), DIMENSION(:,:), ALLOCATABLE :: gFT
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   nGy = SIZE( g, 1 )
   nGx = SIZE( g, 2 )

   xCenter = x1 + dx*REAL( nGx/2, KIND=pr)
   yCenter = y1 + dy*REAL( nGY/2, KIND=pr)

   Ny = SIZE( gInt, 1 )
   Nx = SIZE( gInt, 2 )

   xConst = 1.0_pr/(REAL( nGx, KIND=pr ))
   yConst = 1.0_pr/(REAL( nGy, KIND=pr ))

   xShift = xCenter - ( xFirst + deltaX*REAL( Nx/2, KIND=pr ) )
   yShift = yCenter - ( yFirst + deltaY*REAL( Ny/2, KIND=pr ) )

   dfx = 1.0_pr/(REAL( nGx, KIND=pr )*dx)
   dfy = 1.0_pr/(REAL( nGy, KIND=pr )*dy)

   ALLOCATE( xPhasor(nGx), yPhasor(nGy), gFT(nGy,nGx) )

   xPhasor = xConst * EXP( -TWOPI*i*xShift*CoordinateVector( nGx, dfx ) )
   yPhasor = yConst * EXP( -TWOPI*i*yShift*CoordinateVector( nGy, dfy ) )

   gFT = g
   CALL FFTshift( gFT, inverse=.TRUE. )
   CALL FFT2D( gFT )
   CALL FFTshift( gFT )
   DO k = 1,nGx
   	DO j = 1,nGy
         gFT(j,k) = yPhasor(j)*xPhasor(k)*gFT(j,k)
      END DO
   END DO

   DEALLOCATE( xPhasor, yPhasor )

   CALL CZT2D( gFT, dfy, dfx, deltaY, deltaX, inverse=.TRUE., transform=gInt )

   DEALLOCATE( gFT )

   END SUBROUTINE TrigInterp2D

!DEC$ ENDIF

END MODULE Array_Routines
