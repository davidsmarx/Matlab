!     Last change:  RSB  15 Nov 2000    1:49 pm
MODULE Kinds

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   INTEGER, PARAMETER :: sp  = KIND( 1.0 )            ! Single precision, real
   INTEGER, PARAMETER :: spc = KIND( (1.0, 1.0) )     ! Single precision, complex
   INTEGER, PARAMETER :: dp  = KIND( 1.0d0 )          ! Double precision, real
   INTEGER, PARAMETER :: dpc = KIND( (1.0d0, 1.0d0) ) ! Double precision, complex

   INTEGER, PARAMETER :: pr  = dp      ! Default precision, real
   INTEGER, PARAMETER :: prc = dpc     ! Default precision, complex

END MODULE Kinds
