!     Last change:  RSB   5 Mar 2002    5:10 pm
MODULE Utility_Routines

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

   PUBLIC
   PRIVATE :: WriteToDisk1D, WriteToDisk2D

   INTERFACE WriteToDisk
      MODULE PROCEDURE WriteToDisk1D,  &
                       WriteToDisk2D
   END INTERFACE

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION FormatDate( dateVal ) RESULT( date )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! From the character variable dateVal (as output by DATE_AND_TIME), the
   ! year, month and day are extracted and a date of the form yyyy/mm/dd is
   ! returned.

   !      6 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

	CHARACTER(LEN=8), INTENT(IN) :: dateVal
	CHARACTER(LEN=10)            :: date

	! Local variables
	CHARACTER(LEN=4) :: year
	CHARACTER(LEN=2) :: month
	CHARACTER(LEN=2) :: day
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	year  = dateVal(1:4)
	month = dateVal(5:6)
	day   = dateVal(7:8)

	date = year//"/"//month//"/"//day

	END FUNCTION FormatDate



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION FormatTime( timeVal ) RESULT( time )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! From the character variable timeVal (as output by DATE_AND_TIME), the
   ! hour, minute and second are extracted and a time of the form hh:mm:ss is
   ! returned. Fractions of a second are ignored.

   !      6 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

	CHARACTER(LEN=10), INTENT(IN) :: timeVal
	CHARACTER(LEN=8)              :: time

	! Local variables
	CHARACTER(LEN=2) :: hours
	CHARACTER(LEN=2) :: minutes
	CHARACTER(LEN=2) :: seconds
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	hours   = timeVal(1:2)
	minutes = timeVal(3:4)
	seconds = timeVal(5:6)

	time = hours//":"//minutes//":"//seconds

	END FUNCTION FormatTime



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE FormattedTimeAndDate( time, date )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calls DATE_AND_TIME, and returns time and date in the form
   !
   !	time -> hh:mm:ss
   !
   !	date -> yyyy/mm/dd

   !      6 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

	CHARACTER(LEN=8),  INTENT(OUT) :: time
	CHARACTER(LEN=10), INTENT(OUT) :: date

	INTRINSIC DATE_AND_TIME

	! Local variables
	CHARACTER(LEN=8)  :: dateVal
	CHARACTER(LEN=10) :: timeVal
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   CALL DATE_AND_TIME( dateVal, timeVal )

	time = FormatTime( timeVal )
	date = FormatDate( dateVal )

	END SUBROUTINE FormattedTimeAndDate



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE WriteToDisk2D( array, fileName )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Writes a formatted file with specified name (with ".txt" appended). The
   ! array is written out by rows; the format is large enough for 2048 values.

   !     17 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),        DIMENSION(:,:), INTENT(IN) :: array
   CHARACTER(LEN=*),                INTENT(IN) :: fileName

   INTRINSIC SIZE, TRIM

   ! Local variables
   INTEGER             :: j, k
   INTEGER, PARAMETER  :: unit = 11
   INTEGER             :: ioFlag
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   OPEN ( unit, FILE= TRIM( fileName )//".txt", &
                STATUS="REPLACE",               &
                FORM="FORMATTED",               &
                ACTION="WRITE",                 &
                POSITION="REWIND",              &
                IOSTAT= ioFlag )
   IF ( ioFlag > 0 ) CALL ErrorExit( "WriteToDisk: error opening output unit" )
   DO j = 1,SIZE( array, 1 )
      WRITE ( unit, "(1X,2048ES12.4E2)" )  &
            ( array(j,k), k=1,SIZE( array, 2 ) )
   END DO
   CLOSE ( unit )

   END SUBROUTINE WriteToDisk2D



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE WriteToDisk1D( array, fileName )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Writes a formatted file with specified name (with ".txt" appended). The
   ! format is large enough for 2048 values.

   !     17 Nov  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   REAL(pr),        DIMENSION(:), INTENT(IN) :: array
   CHARACTER(LEN=*),              INTENT(IN) :: fileName

   INTRINSIC SIZE, TRIM

   ! Local variables
   INTEGER             :: j
   INTEGER, PARAMETER  :: unit = 11
   INTEGER             :: ioFlag
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   OPEN ( unit, FILE= TRIM( fileName )//".txt", &
                STATUS="REPLACE",               &
                FORM="FORMATTED",               &
                ACTION="WRITE",                 &
                POSITION="REWIND",              &
                IOSTAT= ioFlag )
   IF ( ioFlag > 0 ) CALL ErrorExit( "WriteToDisk: error opening output unit" )

   WRITE ( unit, "(1X,2048ES12.4E2)" ) ( array(j), j=1,SIZE( array ) )

   CLOSE ( unit )

   END SUBROUTINE WriteToDisk1D



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE WaitForResponse( message, action )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Writes a message to standard output and waits for user response; default,
   ! if action is not present, is "continue"

   !     18 Dec  '00
   !      1 Mar  '02    action input variable is optional
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   CHARACTER(LEN=*), INTENT(IN)           :: message
   CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: action    ! Program stop if "quit"

   INTRINSIC TRIM, PRESENT

   ! Local variable:
   CHARACTER(LEN=1)  :: ignoredCharacter
   CHARACTER(LEN=32) :: what
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   what = "continue"
   IF ( PRESENT( action ) ) what = action

   WRITE ( * , "(1X,A)", ADVANCE = "NO" )            &
          message//"; press <Enter> to "//TRIM( what )
   READ  ( * , "(1X,A)" ) ignoredCharacter

   SELECT CASE ( what )
      CASE ( "quit", "QUIT", "Quit" )
         STOP
      CASE ( "continue", "Continue", "CONTINUE" )
         RETURN
      CASE DEFAULT
   END SELECT

   END SUBROUTINE WaitForResponse

END MODULE Utility_Routines
