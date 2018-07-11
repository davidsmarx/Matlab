!     Last change:  RSB  23 Sep 2002    4:46 pm
MODULE Polygons

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

   PRIVATE

   PUBLIC :: ASSIGNMENT ( = )
   PUBLIC :: AddVertex,        ShiftPolygon,     ClearPolygon
   PUBLIC :: RotatePolygon,    RotatePolygon_3d
   PUBLIC :: PolygonArea,      PolygonOverlap,   PolygonOrientation
   PUBLIC :: MakeHexagon,      MakeRectangle
   PUBLIC :: NumberOfVertices, GetVertex,        PolygonIsConvex
   PUBLIC :: BoundingRectangle

   ! DEFINITION OF POLYGON OBJECT..............................................
   TYPE, PUBLIC :: polygon
      PRIVATE
      INTEGER                           :: numVertices = 0
      REAL(pr), DIMENSION(:,:), POINTER :: vertex => NULL()
   END TYPE polygon
   !...........................................................................

   INTERFACE ASSIGNMENT ( = )
      MODULE PROCEDURE InitializePolygonToArray
      MODULE PROCEDURE SetPolygonEqualToPolygon
   END INTERFACE

   INTEGER, PUBLIC,  PARAMETER :: notConvex =  0
   INTEGER, PUBLIC,  PARAMETER :: clockwise = -1

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE InitializePolygonToArray( poly, vertexArray )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Assign an array of vertex coordinates to a polygon. The array must have
   !  size 2xN, where N >= 3. The form of the array is
   !                 x1 x2 x3 ... <- first row,  x-coordinates
   !                 y1 y2 y3 ... <- second row, y-coordinates
   !     12 Oct  '00
   !      1 Jun  '01    Test for allocated poly%vertex array
   !     19 July '01    Use (:,:) when assigning poly%vertex array
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE( polygon ),          INTENT(OUT) :: poly
   REAL(pr), DIMENSION(:,:), INTENT(IN)  :: vertexArray

   INTRINSIC SIZE, ASSOCIATED

   ! Local variables
   INTEGER n1, n2
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   n1 = SIZE( vertexArray, 1 )
   n2 = SIZE( vertexArray, 2 )

   IF ( n1 /= 2  .OR.  n2 < 3 ) THEN
      CALL ErrorExit( "InitializePolygon: incorrect input for polygon" )
   END IF

   IF ( ASSOCIATED( poly%vertex ) ) DEALLOCATE( poly%vertex )
   poly%numVertices = n2
   ALLOCATE( poly%vertex(2,poly%numVertices) )
   poly%vertex = vertexArray

   END SUBROUTINE InitializePolygonToArray



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SetPolygonEqualToPolygon( poly, poly2 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Set the output polygon, poly, equal to the input, poly2
   !     16 Oct  '00
   !     28 May  '01    Number of vertices always given by poly%numVertices
   !     19 July '01    Use (:,:) when assigning poly%vertex array
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   TYPE( polygon ), INTENT(OUT) :: poly
   TYPE( polygon ), INTENT(IN ) :: poly2

   INTRINSIC ASSOCIATED

   ! Local variables
   INTEGER            :: flag
   INTEGER, PARAMETER :: ok = 0
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ASSOCIATED( poly%vertex ) ) THEN
      IF ( poly%numVertices /= poly2%numVertices ) THEN
         DEALLOCATE( poly%vertex, STAT = flag )
         IF ( flag /= ok ) &
            CALL ErrorExit( "SetPolygonEqualToPolygon: deallocation error" )
         ALLOCATE( poly%vertex(2,poly2%numVertices), STAT = flag )
         IF ( flag /= ok ) &
            CALL ErrorExit( "SetPolygonEqualToPolygon: allocation error 1" )
      END IF
   ELSE
      ALLOCATE( poly%vertex(2,poly2%numVertices), STAT = flag )
      IF ( flag /= ok ) &
         CALL ErrorExit( "SetPolygonEqualToPolygon: allocation error 2" )
   END IF

   poly%numVertices = poly2%numVertices
   poly%vertex(:,:) = poly2%vertex(:,:)

   END SUBROUTINE SetPolygonEqualToPolygon



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ShiftPolygon( poly, xc, yc  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Displaces all vertices by specified distances in x and y

   !     26 Oct  '00
   !     28 May '01     Number of vertices always given by poly%numVertices
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(polygon), INTENT(IN OUT) :: poly
   REAL(pr),      INTENT(IN)     :: xc, yc

   INTRINSIC ABS, EPSILON

   ! Local variables
   INTEGER  :: k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( xc ) > EPSILON(0.0_pr)  .OR.   ABS( yc ) > EPSILON(0.0_pr) ) THEN

      DO k = 1,poly%numVertices
         poly%vertex(1,k) = poly%vertex(1,k) + xc
         poly%vertex(2,k) = poly%vertex(2,k) + yc
      END DO

   END IF

   END SUBROUTINE ShiftPolygon



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE RotatePolygon( poly, angle, xc, yc  )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Rotates the polygon counterclockwise about the point xc, yc
   ! (default is xc = yc = 0 if these parameters are not present)

   !     24 Sept '01
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   TYPE(polygon), INTENT(IN OUT)       :: poly
   REAL(pr),      INTENT(IN)           :: angle
   REAL(pr),      INTENT(IN), OPTIONAL :: xc, yc

   INTRINSIC ABS, EPSILON, PRESENT, MATMUL

   ! Local variables
   REAL(pr) :: xRotCenter, yRotCenter
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( angle ) < EPSILON(0.0_pr) ) RETURN

   xRotCenter = 0.0_pr
   yRotCenter = 0.0_pr
   IF ( PRESENT( xc ) ) xRotCenter = xc
   IF ( PRESENT( yc ) ) yRotCenter = yc

   CALL ShiftPolygon( poly, -xRotCenter, -yRotCenter )
   poly%vertex(:,:) = MATMUL( R_2d( angle ), poly%vertex(:,:) )
   CALL ShiftPolygon( poly, +xRotCenter, +yRotCenter )

   END SUBROUTINE RotatePolygon



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE RotatePolygon_3d( poly, angle, axis )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Rotates the polygon counterclockwise about an axis defined by
   ! the input 3-vector, "axis." The vertices of "poly" are assumed
   ! to initially lie in the plane z = 0. After rotation, the output
   ! polygon is the projection of the rotated polygon onto the plane z = 0.

   !     19 Sept '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   TYPE(polygon), INTENT(IN OUT)           :: poly
   REAL(pr),      INTENT(IN)               :: angle
   REAL(pr),      INTENT(IN), DIMENSION(3) :: axis

   INTRINSIC ABS, EPSILON, MATMUL

   ! Local variables
   INTEGER :: n
   REAL(pr), DIMENSION(3)   :: xyz
   REAL(pr), DIMENSION(3,3) :: R3d
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ABS( angle ) < EPSILON(0.0_pr) ) RETURN

   R3d = General_RotationMatrix( angle, axis )

   DO n = 1,poly%numVertices
      xyz = MATMUL( R3d, (/ poly%vertex(1,n), poly%vertex(2,n), 0.0_pr /) )
      poly%vertex(1,n) = xyz(1)
      poly%vertex(2,n) = xyz(2)
   END DO

   END SUBROUTINE RotatePolygon_3d



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ClearPolygon( poly )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Deallocate storage for vertices
   !     24 Oct  '00
   !      5 Dec  '01    Removed unused local variable
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(polygon), INTENT(IN OUT) :: poly

   INTRINSIC ASSOCIATED

   ! Local variables (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   poly%numVertices = 0
   IF ( ASSOCIATED( poly%vertex ) )  DEALLOCATE( poly%vertex )

   END SUBROUTINE ClearPolygon




  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION PolygonIsConvex( inputPolygon ) RESULT( testResult )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Returns .TRUE. if input is a valid convex polygon; otherwise .FALSE.
   !  13 Oct  '00
   !  17 Apr  '01    Corrected testing of orientation when more than 2 vertices
   !   1 Jun  '01    Test for undefined polygon
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   TYPE( polygon ), INTENT(IN) :: inputPolygon
   LOGICAL                     :: testResult

   INTRINSIC ASSOCIATED

   ! Local variables
   ! (none)
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( .NOT.ASSOCIATED( inputPolygon%vertex ) ) THEN
      CALL ErrorExit( "PolygonIsConvex: input polygon is not defined" )
   END IF

   IF ( inputPolygon%numVertices < 3  .OR.              &
        PolygonOrientation( inputPolygon ) == notConvex )  THEN
      testResult = .FALSE.
   ELSE
      testResult = .TRUE.
   END IF

   END FUNCTION PolygonIsConvex



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE AddVertex( poly, newVertex )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Adds a new vertex to the input polygon (at the "end" of the array);
   !  polygons can be built vertex by vertex.
   !     24 Oct  '00
   !     28 May  '01    Completely revised
   !     19 July '01    Use (:,:) when assigning poly%vertex array
   !      2 Aug  '01    NO  (:,:)  "       "         "         "
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   TYPE( polygon ),        INTENT(IN OUT) :: poly
   REAL(pr), DIMENSION(2), INTENT(IN)     :: newVertex

   INTRINSIC MAX, ASSOCIATED

   !  Local variables
   REAL(pr), DIMENSION(:,:), ALLOCATABLE :: vertexArray
   INTEGER                               :: nold, n
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   nold = poly%numVertices
   poly%numVertices = nold + 1
   n = MAX( nold+1, 2 )                ! Make vertex array at least 2x2

   ALLOCATE( vertexArray(2,n) )
   vertexArray = 0.0_pr

   IF ( .NOT.ASSOCIATED( poly%vertex) ) THEN
      vertexArray(:,1) = newVertex
   ELSE
      vertexArray(:,1:nold) = poly%vertex(:,1:nold)
      vertexArray(:,n)      = newVertex
      DEALLOCATE( poly%vertex )
   END IF

   ALLOCATE( poly%vertex(2,n ) )
   poly%vertex = vertexArray

   DEALLOCATE( vertexArray )

   END SUBROUTINE AddVertex


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION NumberOfVertices( inputPolygon ) RESULT( n )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Return the number of vertices for the input polygon
   !     13 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   TYPE( polygon ), INTENT(IN) :: inputPolygon
   INTEGER                     :: n
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   n = inputPolygon%numVertices

   END FUNCTION NumberOfVertices



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION GetVertex( n, inputPolygon ) RESULT( thisVertex )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Return the vertex coordinates for the input polygon
   !     13 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   INTEGER,         INTENT(IN) :: n
   TYPE( polygon ), INTENT(IN) :: inputPolygon
   REAL(pr), DIMENSION(2)      :: thisVertex
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( n > 0  .AND. n <= inputPolygon%numVertices ) THEN
      thisVertex = inputPolygon%vertex(:,n)
   ELSE
      CALL ErrorExit( "GetVertex: requested vertex number is out of range" )
   END IF

   END FUNCTION GetVertex



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION MakeHexagon( side, angle, xc, yc ) RESULT( hexagon )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Generates a 2x6 matrix of vertex coordinates for a hexagon with center
   !  located at (xc,yc), rotated counterclockwise about its center by "angle"
   !	(in radians).

   !     13 July '99
   !     10 Oct  '00    Fortran 95 version
   !     11 Oct  '00    Uses user-defined type "polygon"
   !      1 Jun  '01    Calls MakeArray to generate 2x6 vertexArray
   !     19 July '01    Use (:,:) when assigning hexagon%vertex array
   !      2 Aug  '01    NO  (:,:)  "       "         "         "
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines

   IMPLICIT NONE

   REAL(pr), INTENT(IN)           :: side
   REAL(pr), INTENT(IN), OPTIONAL :: angle
   REAL(pr), INTENT(IN), OPTIONAL :: xc, yc
   TYPE( polygon )                :: hexagon

   INTRINSIC PRESENT, MATMUL, ASSOCIATED

   !  Local variables
   REAL(pr)                 :: L1, L2, L3
   REAL(pr), DIMENSION(2,6) :: vertexArray
   REAL(pr), DIMENSION(2)   :: center
   INTEGER                  :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   L1 = side
   L2 = 0.5_pr*side
   L3 = SQRT3*L2

   vertexArray = MakeArray( 2, 6, (/ -L1,   -L2,   L2,  L1,     L2, -L2,   &
                                    0.0_pr, -L3,  -L3, 0.0_pr,  L3,  L3 /) )
   IF ( PRESENT( angle ) ) THEN
      vertexArray = MATMUL( R_2d( angle ), vertexArray )
   END IF
   IF ( PRESENT( xc )  .AND.  PRESENT( yc ) ) THEN
      center = (/ xc, yc /)
      DO j = 1,2
         DO k = 1,6
            vertexArray(j,k) = vertexArray(j,k) + center(j)
         END DO
      END DO
   END IF

   IF ( .NOT.ASSOCIATED( hexagon%vertex ) ) ALLOCATE( hexagon%vertex(2,6) )
   hexagon%vertex = vertexArray
   hexagon%numVertices = 6

   END FUNCTION MakeHexagon



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION MakeRectangle( length, width, angle, xc, yc ) RESULT( rect )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Generates a 2x4 matrix of vertex coordinates for a rectangle with center
   !  located at (xc,yc), rotated at angle "angle" (radians).

   !     13 July '99
   !     10 Oct  '00    Fortran 95 version
   !     11 Oct  '00    Uses user-defined type "polygon"
   !      1 Jun  '01    Calls MakeArray to generate 2x4 rectArray
   !     19 July '01    Use (:,:) when assigning rect%vertex array
   !      2 Aug  '01    NO  (:,:)  "       "         "         "
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   REAL(pr), INTENT(IN)           :: length
   REAL(pr), INTENT(IN)           :: width
   REAL(pr), INTENT(IN), OPTIONAL :: angle
   REAL(pr), INTENT(IN), OPTIONAL :: xc, yc
   TYPE( polygon )                :: rect

   INTRINSIC PRESENT, MATMUL, ASSOCIATED

   !  Local variables
   REAL(pr)                 :: L, W
   REAL(pr), DIMENSION(2,4) :: rectArray
   REAL(pr), DIMENSION(2)   :: center
   INTEGER                  :: j, k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   L = 0.5_pr*length
   W = 0.5_pr*width

   rectArray = MakeArray( 2, 4, (/ -L,  L,  L, -L,   &
											  -W, -W,  W,  W /) )
   IF ( PRESENT( angle ) ) THEN
      rectArray = MATMUL( R_2d( angle ), rectArray )
   END IF
   IF ( PRESENT( xc )  .AND.  PRESENT( yc ) ) THEN
      center = (/ xc, yc /)
      DO j = 1,2
         DO k = 1,4
            rectArray(j,k) = rectArray(j,k) + center(j)
         END DO
      END DO
   END IF

   IF ( .NOT.ASSOCIATED( rect%vertex ) ) ALLOCATE( rect%vertex(2,4) )
   rect%vertex      = rectArray
   rect%numVertices = 4

   END FUNCTION MakeRectangle



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION PolygonArea( poly ) RESULT( area )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the area of the input polygon. Area is positive if the polygon
   ! has counter-clockwise orientation, and negative for clockwise.
   ! Ref: Bronshtein & Semendyayev, Handbook of Mathematics,
   !      Van Nostrand Reinhold, N.Y. (1985), p.199

   !      6 July '01
   !     19 July '01    Tests number of vertices
   !      9 Aug  '01    Fixed having both "n" and "N" as variables (The compiler
   !                    didn't catch this!). Changed N to nTot
   !     19 Sept '02    Corrected lower limit in DO-loop
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE( polygon ), INTENT(IN) :: poly
   REAL(pr)                    :: area

   !  Local variables
   REAL(pr) :: sum
   INTEGER  :: n, nTot
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   nTot = poly%numVertices

   IF ( nTot < 3 ) CALL ErrorExit( "PolygonArea: not valid input polygon" )

   sum = (poly%vertex(1,nTot)-poly%vertex(1,1)) *  &
   		(poly%vertex(2,nTot)+poly%vertex(2,1))

   DO n = 2,nTot
      sum = sum +                                    &
      		(poly%vertex(1,n-1)-poly%vertex(1,n)) *  &
      		(poly%vertex(2,n-1)+poly%vertex(2,n))
   END DO

   area = 0.5_pr*sum

   END FUNCTION PolygonArea



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION PolygonOrientation( poly ) RESULT( orientation )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates the orientation (+1 or -1) of the input polygon

   !  poly           A polygon
   !  orientation    Value +1 for counter-clockwise orientation,
   !                       -1 for clockwise,
   !                    and 0 if polygon is not convex

   !     19 Apr '00     Robert Benson
   !     11 Oct '00     Fortran 95 version

   !  Reference:
   !  I.O. Angell & G. Griffith,
   !  "High-resolution Computer Graphics Using FORTRAN 77"
   !  Halsted Press, N.Y. (1987), chapter 5
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE( polygon ), INTENT(IN) :: poly
   INTEGER                     :: orientation

   !  Local variables
   REAL(pr), DIMENSION(:,:), ALLOCATABLE :: v2
   INTEGER                               :: n, vertex, orient
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Fill out an array of vertices, repeating the first two original vertices
   n  = poly%numVertices
   IF ( n < 3 ) CALL ErrorExit( "TestPolygonOrientation: too few vertices" )

   ALLOCATE( v2(2,n+2) )
   v2(:,1:n) = poly%vertex(:,1:n)
   v2(:,n+1) = poly%vertex(:,1)
   v2(:,n+2) = poly%vertex(:,2)

   ! Test each vertex for orientation; convexity requires same orietation
   ! for all vertices

   vertex = 2
   orientation = TestVertex( vertex )

   DO vertex = 3,n+1
      orient = TestVertex( vertex )
      IF ( orient /= orientation ) THEN
         orientation = notConvex
         RETURN
      END IF
   END DO

   DEALLOCATE( v2 )

   CONTAINS

     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      FUNCTION TestVertex( vertex ) RESULT( vertexOrientation )
     !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      INTEGER, INTENT(IN)      :: vertex
      INTEGER                  :: vertexOrientation

      INTRINSIC NINT, SIGN

      !  Local variables
      REAL(pr), DIMENSION(2,2) :: v3
      REAL(pr)                 :: det
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      ! Make a 2x2 matrix from line segments that form this vertex
      v3(:,1) = v2(:,vertex)   - v2(:,vertex-1)
      v3(:,2) = v2(:,vertex+1) - v2(:,vertex)

      ! The orientation is the sign of the determinant of the matrix:
      ! + for counter-clockwise, - for clockwise
      det = v3(1,1)*v3(2,2) - v3(2,1)*v3(1,2)
      vertexOrientation = NINT( SIGN( 1.0_pr, det ) )

      END FUNCTION TestVertex

   END FUNCTION PolygonOrientation



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION PolygonOverlap( polygonA, polygonB ) RESULT( newPolygon )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Finds the polygon formed by the intersection of polygons A and B
   !     19 Apr '00     Robert Benson
   !     11 Oct '00     Fortran 95 version
   !     30 May '01     Clear newFeasiblePolygon at end of loop over polygonB

   ! INPUT
   !  polygonA, polygonB   Convex polygons...must have same orientation

   ! OUTPUT
   !  newPolygon           Intersection of the input polygons. The returned
   !                       polygon has the same orientation as input polygons

   !  Reference:
   !  I.O. Angell & G. Griffith,
   !  "High-resolution Computer Graphics Using FORTRAN 77"
   !  Halsted Press, N.Y. (1987), chapter 5

   !  As used here (and in the reference), a "feasible" polygon is one which
   !  may be the intersection of polygons A and B. This routine begins by using
   !  polygon A as an initial guess, then slicing off pieces of A that lie
   !  outside B.
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   IMPLICIT NONE

   TYPE( polygon ), INTENT(IN) :: polygonA, polygonB
   TYPE( polygon )             :: newPolygon

   INTRINSIC DOT_PRODUCT, ABS, EPSILON

   !  Local variables
   INTEGER                :: orientation
   TYPE( polygon )        :: feasiblePolygon
   TYPE( polygon )        :: newFeasiblePolygon
   REAL(pr), DIMENSION(2) :: pointB1, pointB2, pointF1, pointF2
   REAL(pr), DIMENSION(2) :: newVertex
   INTEGER                :: na, nb, nf
   INTEGER                :: vertex_of_B, vertex_of_F
   REAL(pr), DIMENSION(2) :: diff, coeffs
   REAL(pr)               :: const
   REAL(pr)               :: d1, a1, d2, a2
   INTEGER                :: numNewPoints
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   orientation = PolygonOrientation( polygonB )
   IF ( orientation == 0 ) THEN
      CALL ErrorExit( "PolygonOverlap: Polygon B is not convex" )
   END IF
   IF ( PolygonOrientation( polygonA ) /= orientation ) THEN
      CALL ErrorExit( "PolygonOverlap: Both polygons must have same orientation" )
   END IF

   na = NumberOfVertices( polygonA )
   nb = NumberOfVertices( polygonB )
   nf = nb

   ! Use one of the sides of polygon B to "slice" the feasible polygon.
   ! The new feasible polygon is the remainder. Then the new feasible polygon
   ! becomes the current feasible polygon, and the next side of B is used for
   ! the slicing. Initially, the feasible polygon is taken to be polygon A.

   feasiblePolygon = polygonA

   pointB1 = GetVertex( nb, polygonB )

   DO vertex_of_B = 1,nb                  ! Loop over sides of polygon B

      pointB2 = GetVertex( vertex_of_B, polygonB )

      ! Calculate coefficients for the line defined by points B1 and B2
      diff   = pointB2 - pointB1
      coeffs = (/ -diff(2),  diff(1) /)
      const  = -DOT_PRODUCT( coeffs, pointB1 )

      ! Test each edge of the feasible polygon (defined by line segment F1-F2):
      ! does line segment B1-B2 cross this edge?

      pointF1 = GetVertex( nf, feasiblePolygon )

      ! This point lies on or inside the slicing edge if "distance" d1
      ! is non-negative. d1 is in units of the length of B1-B2.
      d1 = ( DOT_PRODUCT( coeffs, pointF1 ) + const )*orientation
      a1 = ABS( d1 )
      IF ( a1 < EPSILON(a1) ) THEN
         d1 = 0.0_pr
      END IF

      numNewPoints = 0
      DO vertex_of_F = 1,nf               ! Loop over points of feasible polygon

         ! If point F1 lies on or inside the slicing edge,
         ! add it to the new feasible polygon
         IF ( d1 >= 0.0_pr ) THEN
            numNewPoints = numNewPoints + 1
            CALL AddVertex( newFeasiblePolygon, pointF1 )
         END IF

         pointF2 = GetVertex( vertex_of_F, feasiblePolygon )
         d2 = ( DOT_PRODUCT( coeffs, pointF2 ) + const )*orientation
         a2 = abs( d2 )
         IF ( a2 < EPSILON(a2) ) THEN
            d2 = 0.0_pr
         END IF

         ! If pointF1 and pointF2 are on opposite sides of the slicing edge,
         ! find the point of intersection of the line segments and add it to
         ! the new feasible polygon
         IF ( d1 /= 0.0_pr  .AND.  d2 /= 0.0_pr  .AND.  d1*d2 < 0.0_pr ) THEN
            numNewPoints = numNewPoints + 1
            newVertex    = ( a2*pointF1 + a1*pointF2 )/( a1 + a2 )
            CALL AddVertex( newFeasiblePolygon, newVertex )
         END IF

         ! The second point of this edge of the current feasible polygon will be
         ! the first point of the next edge
         d1      = d2
         a1      = a2
         pointF1 = pointF2

      END DO               ! End of loop over points of current feasible polygon

      ! There is no overlap if the new feasible polygon has less than 3 sides
      IF ( numNewPoints < 3 ) THEN
         CALL ErrorExit( "PolygonOverlap: overlap does not exist" )
      END IF

      ! Now slice the new feasible polygon, using the next side of polygon B
      nf = numNewPoints
      feasiblePolygon = newFeasiblePolygon

      ! Move on to the next edge of polygon B
      pointB1 = pointB2
      CALL ClearPolygon( newFeasiblePolygon )

   END DO                  ! End of loop over sides of polygon B

   newPolygon = feasiblePolygon
   CALL ClearPolygon( feasiblePolygon )

   END FUNCTION PolygonOverlap



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   FUNCTION BoundingRectangle( poly ) RESULT( xyBounds )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Returns the x and y values of the rectangle that encloses the input
   !  polygon
   !     22 Apr '02     Robert Benson

   ! INPUT
   !  polygon		A polygon

   ! OUTPUT
   !  xyBounds		A 2x2 matrix whose first row is x-values and second row
   !              is y-values that enclose the polygon

   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   TYPE( polygon ), INTENT(IN) :: poly

   REAL(pr),    DIMENSION(2,2) :: xyBounds

   INTRINSIC MINVAL, MAXVAL

   ! Local variables
   REAL(pr) :: xMin, xMax, yMin, yMax
   ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   xMin = MINVAL( poly%vertex(1,:) )
   xMax = MAXVAL( poly%vertex(1,:) )
   yMin = MINVAL( poly%vertex(2,:) )
   yMax = MAXVAL( poly%vertex(2,:) )

   xyBounds = MakeArray( 2, 2, (/ xMin, xMax,   &
                                  yMin, yMax /) )

   END FUNCTION BoundingRectangle

END MODULE Polygons
