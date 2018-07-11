!     Last change:  RSB  20 Dec 2001    9:59 am
MODULE Constants

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds

   REAL(pr), PARAMETER :: PI     = 3.141592653589793_pr
   REAL(pr), PARAMETER :: TWOPI  = 6.283185307179586_pr
   REAL(pr), PARAMETER :: HALFPI = 1.570796326794896_pr
   REAL(pr), PARAMETER :: SQRT2  = 1.414213562373095_pr
   REAL(pr), PARAMETER :: SQRT3  = 1.732050807568877_pr

   COMPLEX(prc), PARAMETER :: i = ( 0.0_pr, 1.0_pr )

END MODULE Constants
