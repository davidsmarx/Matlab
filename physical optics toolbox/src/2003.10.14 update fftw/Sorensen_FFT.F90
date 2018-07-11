!     Last change:  RSB  18 Dec 2001    5:43 pm
MODULE Sorensen_FFT

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds
   USE Math_Routines, ONLY : PowerOf2
   USE Error_Exit

   IMPLICIT NONE

   PRIVATE
   PUBLIC FFT1D_sorensen, FFT2D_sorensen

   ! Tables of sines, cosines, and bit-reversal quantities
   INTEGER, SAVE :: N = 0
   INTEGER, SAVE :: M
   REAL(pr), DIMENSION(:), ALLOCATABLE, SAVE :: CT1, ST1, CT3, ST3
   INTEGER,  DIMENSION(:), ALLOCATABLE, SAVE :: ITAB

   ! Real and imaginary parts of the input array
   REAL(pr), DIMENSION(:), ALLOCATABLE :: X, Y

   CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FFT1D_sorensen( cpxArray, inverse )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Power-of-2, one-dimensional, split-radix, complex FFT. Calls CTFFTSR,
   ! originally written by H.V. Sorensen (1987), and slightly modified here
   ! for Fortran 90.
   !     27 Nov '00     R.S. Benson    Lockheed Martin Space Systems Co.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:), INTENT(IN OUT) :: cpxArray
   LOGICAL,      OPTIONAL,     INTENT(IN)     :: inverse

   INTRINSIC SIZE, ALLOCATED, CEILING, SQRT, PRESENT, REAL, AIMAG, CMPLX

   ! Local variables:
   LOGICAL :: inv
   INTEGER :: k
   INTEGER :: nStoTrig, nStoTab
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( N /= SIZE( cpxArray ) ) THEN

      ! Allocate enough storage for trig functions and bit-reversal table
      N = SIZE( cpxArray )
      M = PowerOf2( N )

      IF ( 2**M /= N ) CALL ErrorExit( "FFT1D_sorensen: input dimension "//  &
                                       "not a power of 2" )
      IF ( ALLOCATED( CT1 ) ) THEN
         DEALLOCATE( X, Y, CT1, CT3, ST1, ST3, ITAB )
      END IF
      nStoTrig = N/8 - 1
      nStoTab  = CEILING( SQRT( REAL( 2*N ) ) ) + 1
      ALLOCATE( CT1(nStoTrig), CT3(nStoTrig), ST1(nStoTrig), ST3(nStoTrig),  &
                ITAB(nStoTab) )
      CALL TINIT

      ! Allocate storage for real and imaginary parts of input function
      ALLOCATE( X(N), Y(N) )

   END IF

   IF ( PRESENT( inverse ) ) THEN
      inv = inverse
   ELSE
      inv = .FALSE.
   END IF

   X =  REAL( cpxArray )
   Y = AIMAG( cpxArray )

   IF ( inv ) Y = -Y

   CALL CTFFTSR

   IF ( inv ) Y = -Y

   DO k = 1,N
      cpxArray(k) = CMPLX( X(k), Y(k), KIND=prc )
   END DO

   END SUBROUTINE FFT1D_sorensen



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FFT2D_sorensen( cpxArray, inverse )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Power-of-2, two-dimensional, complex FFT. Calls FFT1D based on CTFFTSR,
   ! originally written by H.V. Sorensen (1987), and slightly modified here
   ! for Fortran 90.
   !     28 Nov '00     R.S. Benson    Lockheed Martin Space Systems Co.
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   COMPLEX(prc), DIMENSION(:,:), INTENT(IN OUT), TARGET :: cpxArray
   LOGICAL,      OPTIONAL,       INTENT(IN)             :: inverse

   INTRINSIC PRESENT, SIZE

   ! Local variables:
   COMPLEX(prc), DIMENSION(:), POINTER :: data
   LOGICAL :: inv
   INTEGER :: k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( PRESENT( inverse ) ) THEN
      inv = inverse
   ELSE
      inv = .FALSE.
   END IF

   DO k = 1, SIZE( cpxArray, 1 )
      data => cpxArray(k,:)
      CALL FFT1D_sorensen( data, inv )
   END DO
   DO k = 1, SIZE( cpxArray, 2 )
      data => cpxArray(:,k)
      CALL FFT1D_sorensen( data, inv )
   END DO
   NULLIFY( data )

   END SUBROUTINE FFT2D_sorensen



!!=================================================================!!
!!                                                                 !!
!!  Subroutine CTFFTSR(X,Y,M,CT1,CT3,ST1,ST3,ITAB):                !!
!!      An in-place, split-radix complex FFT program               !!
!!      Decimation-in-frequency, cos/sin in third loop             !!
!!      and is looked-up in table. Tables CT1,CT3,ST1,ST3          !!
!!      have to have length>=N/8-1. The bit reverser uses partly   !!
!!      table lookup.                                              !!
!!                                                                 !!
!!  Input/output                                                   !!
!!      X    Array of real part of input/output (length >= N)      !!
!!      Y    Array of imaginary part of input/output (length >= N) !!
!!      M    Transform length is N=2**M                            !!
!!      CT1  Array of cos() table (length >= N/8-1)                !!
!!      CT3  Array of cos() table (length >= N/8-1)                !!
!!      ST1  Array of sin() table (length >= N/8-1)                !!
!!      ST3  Array of sin() table (length >= N/8-1)                !!
!!      ITAB Array of bitreversal indices (length >= sqrt(2*N)     !!
!!                                                                 !!
!!  Calls:                                                         !!
!!      CTSTAG                                                     !!
!!      and TINIT has to be called before this program!!           !!
!!                                                                 !!
!!  Author:                                                        !!
!!      H.V. Sorensen,   University of Pennsylvania,  Dec. 1984    !!
!!                       Arpa address: hvs@ee.upenn.edu            !!
!!  Modified:                                                      !!
!!      H.V. Sorensen,   University of Pennsylvania,  Jul. 1987    !!
!!                                                                 !!
!!  Reference:                                                     !!
!!      Sorensen, Heideman, Burrus :"On computing the split-radix  !!
!!      FFT", IEEE Tran. ASSP, Vol. ASSP-34, No. 1, pp. 152-156    !!
!!      Feb. 1986                                                  !!
!!      Mitra&Kaiser: "Digital Signal Processing Handbook, Chap.   !!
!!      8, page 491-610, John Wiley&Sons, 1993                     !!
!!                                                                 !!
!!      This program may be used and distributed freely as long    !!
!!      as this header is included                                 !!
!!                                                                 !!
!!=================================================================!!
      SUBROUTINE CTFFTSR

     ! Local variables:
     INTEGER  :: ITS, N2, N4, IS, ID, M2, NBIT
     INTEGER  :: I, J, K, L, I1, J0
     REAL(pr) :: T1

!-------L shaped butterflies----------------------------------------
      ITS = 1
      N2 = 2*N
      DO K = 1, M-1
         N2 = N2/2
         N4 = N2/4
         CALL CTSTAG(N2,N4,ITS,X(1:N),X(N4+1:N),X(2*N4+1:N),X(3*N4+1:N), &
                               Y(1:N),Y(N4+1:N),Y(2*N4+1:N),Y(3*N4+1:N) )
         ITS = 2 * ITS
      END DO
!-------Length two butterflies--------------------------------------
      IS = 1
      ID = 4
      DO
         DO I1 = IS,N,ID
            T1      = X(I1)
            X(I1)   = T1 + X(I1+1)
            X(I1+1) = T1 - X(I1+1)
            T1      = Y(I1)
            Y(I1)   = T1 + Y(I1+1)
            Y(I1+1) = T1 - Y(I1+1)
         END DO
         IS = 2*ID - 1
         ID = 4*ID
         IF ( IS >= N ) EXIT
      END DO
!-------Digit reverse counter---------------------------------------
      M2 = M/2
      NBIT = 2**M2
      DO K = 2, NBIT
         J0 = NBIT * ITAB(K) + 1
         I = K
         J = J0
         DO L = 2, ITAB(K)+1
            T1   = X(I)
            X(I) = X(J)
            X(J) = T1
            T1   = Y(I)
            Y(I) = Y(J)
            Y(J) = T1
            I = I + NBIT
            J = J0 + ITAB(L)
         END DO
      END DO

      END SUBROUTINE CTFFTSR



!===================================================================!
!   Subroutine CTSTAG - the work-horse of CTFFTSR                   !
!       Computes a stage of a length N split-radix transform        !
!   Author:                                                         !
!       H.V. Sorensen,   University of Pennsylvania,  Jul. 1987     !
!===================================================================!
        SUBROUTINE CTSTAG( N2,N4,ITS,X1,X2,X3,X4,Y1,Y2,Y3,Y4 )

        INTEGER,                INTENT(IN)     :: N2, N4, ITS
        REAL(pr), DIMENSION(:), INTENT(IN OUT) :: X1 ,X2, X3, X4, &
                                                  Y1, Y2, Y3, Y4

       ! Local variables:
        REAL(pr), PARAMETER :: cos45 = 0.7071067818654752440D+0
        INTEGER  :: N8, IS, ID, IT
        INTEGER  :: I1, I2
        INTEGER  :: I, J, JN
        REAL(pr) :: T1, T2, T3, T4, T5

        N8 = N4/2
!-------Zero butterfly----------------------------------------------
        IS = 0
        ID = 2*N2
        DO
           DO I1 = IS+1,N,ID
              T1     = X1(I1) - X3(I1)
              X1(I1) = X1(I1) + X3(I1)
              T2     = Y2(I1) - Y4(I1)
              Y2(I1) = Y2(I1) + Y4(I1)
              X3(I1) = T1     + T2
              T2     = T1     - T2
              T1     = X2(I1) - X4(I1)
              X2(I1) = X2(I1) + X4(I1)
              X4(I1) = T2
              T2     = Y1(I1) - Y3(I1)
              Y1(I1) = Y1(I1) + Y3(I1)
              Y3(I1) = T2     - T1
              Y4(I1) = T2     + T1
           END DO
           IS = 2*ID - N2
           ID = 4*ID
           IF (IS >= N) EXIT
        END DO

        IF  ( N4 <= 1 ) RETURN
!-------N/8 butterfly-----------------------------------------------
        IS = 0
        ID = 2*N2
        DO
           DO I1 = IS+1+N8,N,ID
              T1     = X1(I1) - X3(I1)
              X1(I1) = X1(I1) + X3(I1)
              T2     = X2(I1) - X4(I1)
              X2(I1) = X2(I1) + X4(I1)
              T3     = Y1(I1) - Y3(I1)
              Y1(I1) = Y1(I1) + Y3(I1)
              T4     = Y2(I1) - Y4(I1)
              Y2(I1) = Y2(I1) + Y4(I1)
              T5     = (T4 - T1)*cos45
              T1     = (T4 + T1)*cos45
              T4     = (T3 - T2)*cos45
              T2     = (T3 + T2)*cos45
              X3(I1) = T4 + T1
              Y3(I1) = T4 - T1
              X4(I1) = T5 + T2
              Y4(I1) = T5 - T2
           END DO
           IS = 2*ID - N2
           ID = 4*ID
           IF (IS >= N-1) EXIT
        END DO

        IF  (N8 <= 1) RETURN
!-------General butterfly. Two at a time----------------------------
        IS = 1
        ID = N2*2
        DO
           DO  I = IS, N, ID
              IT = 0
              JN = I + N4
              DO J=1, N8-1
                 IT = IT+ITS
                 I1 = I+J
                 T1     = X1(I1) - X3(I1)
                 X1(I1) = X1(I1) + X3(I1)
                 T2     = X2(I1) - X4(I1)
                 X2(I1) = X2(I1) + X4(I1)
                 T3     = Y1(I1) - Y3(I1)
                 Y1(I1) = Y1(I1) + Y3(I1)
                 T4     = Y2(I1) - Y4(I1)
                 Y2(I1) = Y2(I1) + Y4(I1)
                 T5 = T1 - T4
                 T1 = T1 + T4
                 T4 = T2 - T3
                 T2 = T2 + T3
                 X3(I1) =  T1*CT1(IT) - T4*ST1(IT)
                 Y3(I1) = -T4*CT1(IT) - T1*ST1(IT)
                 X4(I1) =  T5*CT3(IT) + T2*ST3(IT)
                 Y4(I1) =  T2*CT3(IT) - T5*ST3(IT)
                 I2 = JN - J
                 T1     = X1(I2) - X3(I2)
                 X1(I2) = X1(I2) + X3(I2)
                 T2     = X2(I2) - X4(I2)
                 X2(I2) = X2(I2) + X4(I2)
                 T3     = Y1(I2) - Y3(I2)
                 Y1(I2) = Y1(I2) + Y3(I2)
                 T4     = Y2(I2) - Y4(I2)
                 Y2(I2) = Y2(I2) + Y4(I2)
                 T5 = T1 - T4
                 T1 = T1 + T4
                 T4 = T2 - T3
                 T2 = T2 + T3
                 X3(I2) =  T1*ST1(IT) - T4*CT1(IT)
                 Y3(I2) = -T4*ST1(IT) - T1*CT1(IT)
                 X4(I2) = -T5*ST3(IT) - T2*CT3(IT)
                 Y4(I2) = -T2*ST3(IT) + T5*CT3(IT)
              END DO
           END DO
           IS = 2*ID - N2 +1
           ID = 4*ID
           IF (IS >= N) EXIT
        END DO

        END SUBROUTINE CTSTAG



!!=================================================================!!
!!  Subroutine TINIT:                                              !!
!!      Initialize SIN/COS and bit reversal tables                 !!
!!  Author:                                                        !!
!!      H.V. Sorensen,   University of Pennsylvania,  Jul. 1987    !!
!!=================================================================!!
        SUBROUTINE TINIT

        INTRINSIC COS, SIN

       ! Local variables:
       INTEGER  :: I, M2, NBIT, IMAX, LBSS
       REAL(pr) :: ANG, theta, theta3

!-------Sin/Cos table-----------------------------------------------
!       N = 2**M
        ANG = 6.283185307179586D+0/REAL( N, KIND=pr )
        DO I=1,N/8-1
           theta  = ANG * REAL(  I,  KIND=pr )
           theta3 = ANG * REAL( 3*I, KIND=pr )
           CT1(I) = COS(theta)
           CT3(I) = COS(theta3)
           ST1(I) = SIN(theta)
           ST3(I) = SIN(theta3)
        END DO
!-------Bit reversal table------------------------------------------
        M2 = M/2
        NBIT = 2**M2
        IF (2*M2 /= M) M2 = M2 + 1
        ITAB(1) = 0
        ITAB(2) = 1
        IMAX = 1
        DO LBSS = 2, M2
           IMAX = 2 * IMAX
           DO I = 1, IMAX
              ITAB(I)      = 2 * ITAB(I)
              ITAB(I+IMAX) = 1 + ITAB(I)
           END DO
        END DO

        END SUBROUTINE TINIT

END MODULE Sorensen_FFT
