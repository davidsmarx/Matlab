!     Last change:  RSB   7 Jan 2002    9:18 am
MODULE Discrete_Transforms

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds

   IMPLICIT NONE

   PUBLIC
   PRIVATE :: FFTshift2D, FFTshift1D

   INTERFACE FFTshift
      MODULE PROCEDURE FFTshift2D, FFTshift1D
   END INTERFACE

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE CZT2D( array, rowSpacing,     colSpacing,              &
                            rowFreqSpacing, colFreqSpacing, inverse, &
                     transform )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Chirp z-transform of input array.
   !  Ref:  Rabiner, R.W. Schafer & C.M. Rader, "The Chirp z-Transform,"
   !        IEEE Transactions on Audio and Electroacoustics, AU-17, 86-92 (1969)
   !  A one-dimensional version of this transform is described by S. Roose,
   !  B. Brichau & E.W. Stijns, Optics Communications vol. 97, pp.312-318 (1993).
   !  For more detail specific to this implementation, see Tech Note (16 June '00).

   !  Zero frequency is returned at sample point
   !              [ floor(numRowsOut/2)+1, floor(numColsOut/2)+1 ]

   !  The size of the output array is equal to the size of outputArray, if
   !  present; otherwise it is equal to the size of the input array.

   !  Two-dimensional transform done via successive 1-d transforms. Calls an
   !  FFT that does NOT do a 1/N normalization.

   !     17 July '00    R.S. Benson
   !     10 Nov  '00    Fortran 95 version
   !     13 Nov  '00    Calls FFT, rather than fftn
   !     14 Nov  '00    Uses "un-normalized" FFT
   !     27 Nov  '00    Calls revised PowerOf2
   !     31 Jan  '01    Pad to permit use of power-of-2 FFT
   !     17 Dec  '01    Local_Z_Chirp and Local_Z_Power are subroutines

   !  array          Function to be transformed (generally a 2-d array)
   !  rowSpacing     Sample spacing for array, along rows
   !  colSpacing       "       "     "    "      "   columns
   !  rowFreqSpacing Frequency spacing for output, along rows
   !  colFreqSpacing    "         "     "    "       "   coumns
   !  inverse        OPTIONAL input, .TRUE. if inverse transform is desired;
   !                 default is .FALSE.
   !  transform      OPTIONAL output array, required when number of output samples
   !                 is not equal to number of input samples in each direction,
   !                 and also may be specified when the input array should not be
   !                 over-written.
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Math_Routines, ONLY: PowerOf2
   USE Constants
   USE Error_Exit

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT)         :: array
   REAL(pr),                     INTENT(IN)             :: rowSpacing
   REAL(pr),                     INTENT(IN)             :: colSpacing
   REAL(pr),                     INTENT(IN)             :: rowFreqSpacing
   REAL(pr),                     INTENT(IN)             :: colFreqSpacing
   LOGICAL,                      INTENT(IN),     OPTIONAL :: inverse
   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT), OPTIONAL :: transform

   INTRINSIC LOG, MAX, TRANSPOSE, SIZE, PRESENT

   ! Local variables
   LOGICAL :: inv
   INTEGER :: numRowsIn,  numColsIn
   INTEGER :: numRowsOut, numColsOut
   INTEGER :: Lrows,      Lcols
   REAL(pr):: alphaRows,  alphaCols
   INTEGER :: n1, k1
   INTEGER :: j,  k     ! Loop indices
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE :: ZA
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE :: ZB
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE :: Z_CHIRP
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE :: fft_Z_CHIRP
   INTEGER                                 :: maxRows
   INTEGER                                 :: maxCols
   LOGICAL                                 :: inPlace
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE :: dataIn, dataOut
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   numRowsIn = SIZE( array, 1 )
   numColsIn = SIZE( array, 2 )

   IF ( PRESENT( inverse ) ) THEN
      inv = inverse
   ELSE
      inv = .FALSE.
   END IF

   IF ( PRESENT( transform ) ) THEN
      inPlace    = .FALSE.
      numRowsOut = SIZE( transform, 1 )
      numColsOut = SIZE( transform, 2 )
   ELSE
      inPlace    = .TRUE.
      numRowsOut = numRowsIn
      numColsOut = numColsIn
   END IF

   ! Pad to avoid aliasing during convolution.
   ! Use enough padding to permit using a power-of-2 FFT.
   Lrows = 2**PowerOf2( numRowsIn + numRowsOut )
   Lcols = 2**PowerOf2( numColsIn + numColsOut )

   IF ( Lrows <= 0  .OR.  Lcols <= 0 )  &
      CALL ErrorExit( "Lrows or Lcols error in CZT2D" )

   alphaRows = PI*rowFreqSpacing*rowSpacing;
   alphaCols = PI*colFreqSpacing*colSpacing;

   ! Do transform for each column of input array

   n1 = numRowsIn/2
   k1 = numRowsOut/2

   maxRows = MAX( numRowsIn, numRowsOut )
   maxCols = MAX( numColsIn, numColsOut )

   ALLOCATE( ZA(numRowsIn) )
   ALLOCATE( ZB(numRowsOut) )
   ALLOCATE( Z_CHIRP(maxRows) )
   ALLOCATE( fft_Z_CHIRP(Lrows) )
   ALLOCATE( dataIn(numRowsIn), dataOut(numRowsOut) )

   CALL Local_Z_Power( alphaRows, numRowsIn,  2*k1, n1*k1, ZA )
   CALL Local_Z_Power( alphaRows, numRowsOut, 2*n1, n1*k1, ZB )
   CALL Local_Z_Chirp( alphaRows, maxRows, Lrows, Z_CHIRP )

   fft_Z_CHIRP = (/ ( ( 0.0_pr, 0.0_pr ), k = 1,Lrows ) /)

   DO k = 1,numRowsOut
      fft_Z_CHIRP(k) = Z_CHIRP(k)
   END DO
   DO k = Lrows-numRowsIn+2,Lrows
      fft_Z_CHIRP(k) = Z_CHIRP(Lrows+2-k)
   END DO

   CALL FFT1D( fft_Z_CHIRP )

   IF ( inPlace ) THEN
      DO k = 1,numColsIn
         dataIn(:) = array(:,k)
         CALL Local_CZT_1D( dataIn, numRowsIn, numRowsOut, Lrows, &
                              ZA, ZB, fft_Z_CHIRP, inv, dataOut )
         array(:,k) = dataOut(:)
      END DO
   ELSE
      DO k = 1,numColsIn
         dataIn(:) = array(:,k)
         CALL Local_CZT_1D( dataIn, numRowsIn, numRowsOut, Lrows, &
                              ZA, ZB, fft_Z_CHIRP, inv, dataOut )
         transform(:,k) = dataOut(:)
      END DO
   END IF

   ! Do transforms for each row

   n1 = numColsIn/2
   k1 = numColsOut/2

   IF ( numColsIn /= numRowsIn ) THEN
      DEALLOCATE( ZA, dataIn )
      ALLOCATE( ZA(numColsIn), dataIn(numColsIn) )
   END IF
   IF ( numColsOut /= numRowsOut ) THEN
      DEALLOCATE( ZB, dataOut )
      ALLOCATE( ZB(numColsOut), dataOut(numColsOut) )
   END IF
   IF ( maxCols /= maxRows ) THEN
      DEALLOCATE( Z_CHIRP )
      ALLOCATE( Z_CHIRP(maxCols) )
   END IF
   IF ( Lcols /= Lrows ) THEN
      DEALLOCATE( fft_Z_CHIRP )
      ALLOCATE( fft_Z_CHIRP(Lcols) )
   END IF

   CALL Local_Z_Power( alphaCols, numColsIn,  2*k1, n1*k1, ZA )
   CALL Local_Z_Power( alphaCols, numColsOut, 2*n1, n1*k1, ZB )
   CALL Local_Z_Chirp( alphaCols, maxCols, Lcols, Z_CHIRP )

   fft_Z_CHIRP = (/ ( ( 0.0_pr, 0.0_pr ), k = 1,Lcols ) /)

   DO k = 1,numColsOut
      fft_Z_CHIRP(k) = Z_CHIRP(k)
   END DO
   DO k = Lcols-numColsIn+2,Lcols
      fft_Z_CHIRP(k) = Z_CHIRP(Lcols+2-k)
   END DO

   DEALLOCATE( Z_CHIRP )

   CALL FFT1D( fft_Z_CHIRP )

   IF ( inPlace ) THEN
      DO j = 1,numRowsOut
         dataIn(:) = array(j,:)
         CALL Local_CZT_1D( dataIn, numColsIn, numColsOut, Lcols, &
                              ZA, ZB, fft_Z_CHIRP, inv, dataOut )
         array(j,:) = dataOut(:)
      END DO
   ELSE
      DO j = 1,numRowsOut
         dataIn(:) = transform(j,:)
         CALL Local_CZT_1D( dataIn, numColsIn, numColsOut, Lcols, &
                              ZA, ZB, fft_Z_CHIRP, inv, dataOut )
         transform(j,:) = dataOut(:)
      END DO
   END IF

   DEALLOCATE( ZA, ZB, fft_Z_CHIRP )
   DEALLOCATE( dataIn, dataOut )

   CONTAINS


     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE Local_Z_Power( alpha, n, two_k_half, nk1, z_phasor )
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      REAL(pr),                   INTENT(IN)    :: alpha
      INTEGER,                    INTENT(IN)    :: n, two_k_half, nk1
      COMPLEX(prc), DIMENSION(n), INTENT(INOUT) :: z_phasor

      INTRINSIC EXP

     ! Local variables
      INTEGER :: k
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      DO k = 1,n
         z_phasor(k) =  EXP( -i*alpha*( (k-two_k_half-1 )*( k-1 ) + nk1 ) )
      END DO

      END SUBROUTINE Local_Z_Power



     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE Local_Z_Chirp( alpha, n, L, z_chirp )
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      IMPLICIT NONE

      REAL(pr),                   INTENT(IN)    :: alpha
      INTEGER,                    INTENT(IN)    :: n, L
      COMPLEX(prc), DIMENSION(n), INTENT(INOUT) :: z_chirp

      INTRINSIC EXP

     ! Local variables
      INTEGER  :: k
      REAL(pr) :: factor
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Factor required because FFT1D (or its inverse) is not normalized
      factor = (1.0_pr/REAL( L, KIND=pr ))

      DO k = 1,n
         z_chirp(k) = factor * EXP( i*alpha*(k-1)**2 )
      END DO

      END SUBROUTINE Local_Z_CHIRP



     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      SUBROUTINE Local_CZT_1D( array, numIn, numOut, L, za, zb, fft_z_chirp, &
                               invCZT, czt )
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Performs 1-d chirp-z transform on input array. Assumes vectors za, zb
      ! and fft_z_chirp all have the proper sizes. Calls the Singleton FFT.

      INTEGER,                         INTENT(IN)  :: numIn, numOut, L
      COMPLEX(prc), DIMENSION(numIn),  INTENT(IN ) :: array
      COMPLEX(prc), DIMENSION(numIn),  INTENT(IN)  :: za
      COMPLEX(prc), DIMENSION(numOut), INTENT(IN)  :: zb
      COMPLEX(prc), DIMENSION(L),      INTENT(IN)  :: fft_z_chirp
      LOGICAL,                         INTENT(IN)  :: invCZT
      COMPLEX(prc), DIMENSION(numOut), INTENT(OUT) :: czt

      INTRINSIC SHAPE, CONJG

      ! Local variables
      COMPLEX(prc), DIMENSION(L) :: u
      INTEGER :: k
     ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      u = (/ ( ( 0.0_pr, 0.0_pr ), k = 1,L ) /)

      IF ( invCZT ) THEN
         DO k = 1,numIn
            u(k) = CONJG( array(k) ) * za(k)
         END DO
      ELSE
         DO k = 1,numIn
            u(k) = array(k) * za(k)
         END DO
      END IF

      ! Do discrete convolution via FFTs

      CALL FFT1D( u )
      u = u * fft_z_chirp
      CALL FFT1D( u, inverse = .TRUE. )

      DO k = 1,numOut
         czt(k) = u(k) * zb(k)
      END DO

      IF ( invCZT ) czt = CONJG( czt )

      END SUBROUTINE Local_CZT_1D

   END SUBROUTINE CZT2D



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FFT2D( a, inverse )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calls for a 2-d FFT of the input array.

   !      3 Nov  '00
   !     29 Nov  '00    FFT alternatives are Singleton (arbitrary N),
   !                                      or Sorensen (power of 2)
   !     31 Jan  '01    Tests for power-of-2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Sorensen_FFT
   USE Singleton_FFT
   USE Math_Routines, ONLY: PowerOf2

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: a
   LOGICAL,      OPTIONAL,       INTENT(IN)     :: inverse

   INTRINSIC PRESENT, SIZE, SHAPE

   ! Local variables
   LOGICAL :: doInverse
   INTEGER :: Nrows, Ncols
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( PRESENT( inverse ) ) THEN
      doInverse = inverse
   ELSE
      doInverse = .FALSE.
   END IF

   Nrows = SIZE( a, 1 )
   Ncols = SIZE( a, 2 )

   IF ( Nrows == 2**PowerOf2( Nrows )  .AND.  &
        Ncols == 2**PowerOf2( Ncols ) ) THEN
      CALL FFT2D_sorensen( a, doInverse )
   ELSE
      ! The following calls the Singleton FFT
      CALL fftn( a, SHAPE( a ), inv = doInverse )
   END IF

   END SUBROUTINE FFT2D



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FFT1D( a, inverse )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calls for a 1-d FFT of the input array.

   !     29 Nov  '00    FFT alternatives are Singleton (arbitrary N),
   !                                      or Sorensen (power of 2)
   !     31 Jan  '01    Tests for power-of-2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Sorensen_FFT
   USE Singleton_FFT
   USE Math_Routines, ONLY: PowerOf2

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:),INTENT(IN OUT) :: a
   LOGICAL,      OPTIONAL,    INTENT(IN)     :: inverse

   INTRINSIC PRESENT, SIZE, SHAPE

   ! Local variables
   LOGICAL :: doInverse
   INTEGER :: N
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( PRESENT( inverse ) ) THEN
      doInverse = inverse
   ELSE
      doInverse = .FALSE.
   END IF

   N = SIZE( a )

   IF ( N == 2**PowerOf2( N ) ) THEN
      CALL FFT1D_sorensen( a, doInverse )
   ELSE
      ! The following calls the Singleton FFT
      CALL fftn( a, SHAPE( a ), inv = doInverse )
   END IF

   END SUBROUTINE FFT1D



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FFTshift1D( vector, inverse )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Swaps upper and lower halves of 1-d array. If "inverse" is .FALSE.
   ! (the default condition) the values at k = 1,2,...,N-(N/2) are moved to
   ! k = N-(N/2),...,N, while k = N-(N/2)+1,...,N are moved to k = 1,...,(N/2).
   ! For inverse = .FALSE. the shift is in the opposite direction.

   !      3 Nov  '00
   !      4 Jan  '02    Rewritten to correctly treat both even and odd N
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:), INTENT(IN OUT) :: vector
   LOGICAL,      OPTIONAL,     INTENT(IN)     :: inverse

   INTRINSIC SIZE, PRESENT, CSHIFT

   ! Local variables
   INTEGER :: Nshift
   LOGICAL :: inv
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Nshift = (1+SIZE( vector ))/2

   inv = .FALSE.
   IF ( PRESENT( inverse ) ) inv = inverse
   IF ( inv ) Nshift = -Nshift

   vector = CSHIFT( vector, Nshift )

   END SUBROUTINE FFTshift1D



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FFTshift2D( array, inverse )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Two-dimensional version of FFTshift1D. The input array may be large, so
   ! one-dimensional shifts are used in order to minimize requirements for
   ! temporary storage.

   !      8 Nov  '00
   !      7 Jan  '02    Completely rewritten to correctly treat odd and even
   !                    number of samples; calls 1-d CSHIFT to save storage
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT) :: array
   LOGICAL,      OPTIONAL,       INTENT(IN)     :: inverse

   INTRINSIC PRESENT, SIZE, CSHIFT

   ! Local variables
   INTEGER :: Ny, yShift, j
   INTEGER :: Nx, xShift, k
   LOGICAL :: inv
   COMPLEX(prc), DIMENSION(:), ALLOCATABLE :: temp
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Ny = SIZE( array, 1 )
   Nx = SIZE( array, 2 )
   yShift = (Ny+1)/2
   xShift = (Nx+1)/2

   inv = .FALSE.
   IF ( PRESENT( inverse ) ) inv = inverse
   IF ( inv ) THEN
      yShift = -yShift
      xShift = -xShift
   END IF

   ALLOCATE( temp(Ny) )

   ! Shift in y-direction for all x-values
   DO k = 1, Nx
      DO j = 1, Ny
         temp(j) = array(j,k)
      END DO
      temp = CSHIFT( temp, yShift )
      DO j = 1, Ny
         array(j,k) = temp(j)
      END DO
   END DO

   IF ( Nx /= Ny ) THEN
      DEALLOCATE( temp )
      ALLOCATE( temp(Nx) )
   END IF

   ! Shift in x-direction for all y-values
   DO j = 1, Ny
      DO k = 1, Nx
         temp(k) = array(j,k)
      END DO
      temp = CSHIFT( temp, xShift )
      DO k = 1, Nx
         array(j,k) = temp(k)
      END DO
   END DO

   DEALLOCATE( temp )

   END SUBROUTINE FFTshift2D

END MODULE Discrete_Transforms
