!     Last change:  RSB  12 Mar 2001   10:22 am
MODULE Error_Exit

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   IMPLICIT NONE

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ErrorExit( message, action )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Writes error message, then either quits, or waits for user response.

   !     16 Nov  '00
   !      2 Mar  '01    Does not wait for user response
   !     12 Mar  '01    Optional input, "action" (either "wait" or "quit")
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   CHARACTER( LEN = * ), INTENT(IN)           :: message
   CHARACTER( LEN = * ), INTENT(IN), OPTIONAL :: action

   INTRINSIC PRESENT

   ! Local variables
   INTEGER,              PARAMETER :: errUnit = 1
   CHARACTER( LEN = * ), PARAMETER :: errFileName = "C:\Diffraction_Error.txt"
   CHARACTER :: ignoredCharacter
   CHARACTER( LEN = 4 ) :: act
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   WRITE ( * , "(1X,A)" ) message

   IF( PRESENT( action ) ) THEN
      act = action
   ELSE
      act = "quit"
   END IF

   SELECT CASE ( act )
      CASE ( "wait", "WAIT", "Wait" )
         WRITE ( * , "(1X,A)", ADVANCE = "NO" ) "Press <Enter> to quit."
         READ  ( * , "(1X,A)" ) ignoredCharacter
      CASE DEFAULT
         OPEN ( errUnit, FILE   = errFileName, &
                         FORM   = "FORMATTED", &
                         STATUS = "REPLACE" )
         WRITE ( errUnit, "(1X,A)" ) message
         CLOSE ( errUnit )
   END SELECT

   STOP

   END SUBROUTINE ErrorExit


END MODULE Error_Exit
