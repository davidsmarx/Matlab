!     Last change:  RSB  20 Dec 2001   10:57 am
MODULE SI_Units

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds

   REAL(pr), PARAMETER :: meters = 1.0_pr
   REAL(pr), PARAMETER :: centimeters = 1.0e-2_pr * meters
   REAL(pr), PARAMETER :: millimeters = 1.0e-3_pr * meters
   REAL(pr), PARAMETER :: micrometers = 1.0e-6_pr * meters
   REAL(pr), PARAMETER :: nanometers  = 1.0e-9_pr * meters
   REAL(pr), PARAMETER :: picometers  = 1.0e-12_pr * meters
   REAL(pr), PARAMETER :: kilometers  = 1.0e+3_pr * meters
   REAL(pr), PARAMETER :: m  = meters
   REAL(pr), PARAMETER :: cm = centimeters
   REAL(pr), PARAMETER :: mm = millimeters
   REAL(pr), PARAMETER :: nm = nanometers
   REAL(pr), PARAMETER :: pm = picometers
   REAL(pr), PARAMETER :: km = kilometers

END MODULE SI_Units
