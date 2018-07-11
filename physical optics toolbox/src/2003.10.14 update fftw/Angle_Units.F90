!     Last change:  RSB  20 Dec 2001   10:57 am
MODULE Angle_Units

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds
   USE Constants

   REAL(pr), PARAMETER :: radians = 1.0_pr
   REAL(pr), PARAMETER :: milliradians = 1.0e-3_pr
   REAL(pr), PARAMETER :: microradians = 1.0e-6_pr
   REAL(pr), PARAMETER :: nanoradians  = 1.0e-9_pr
   REAL(pr), PARAMETER :: picoradians  = 1.0e-12_pr
   REAL(pr), PARAMETER :: deg = PI/180.0_pr
   REAL(pr), PARAMETER :: arcmin = deg/60.0_pr
   REAL(pr), PARAMETER :: arcsec = arcmin/60.0_pr
   REAL(pr), PARAMETER :: milliarcsec = 1.0e-3_pr*arcsec
   REAL(pr), PARAMETER :: microarcsec = 1.0e-6_pr*arcsec

END MODULE Angle_Units
