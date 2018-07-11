!     Last change:  RSB   1 May 2003    5:11 pm
MODULE Corner_Cubes

!============================
!  Robert S. Benson
!  Lockheed Martin ATC
!  3251 Hanover Street
!  Palo Alto CA 94304
!  robert.s.benson@lmco.com
!============================

   USE Kinds
   USE Polygons
   USE Error_Exit

   PRIVATE
   PUBLIC :: NewCornerCube
   PUBLIC :: SetTilt
   PUBLIC :: SetClockAngle
   PUBLIC :: SetDihedralErrors
   PUBLIC :: ShiftCornerCubeVertex
   PUBLIC :: CornerCubeReflect

   ! DEFINITION OF CORNER CUBE OBJECT..........................................
   ! Components initialized, in accord with Fortran 95
   ! 15 Apr '03   Each gap has separate width value
   !  1 May '03   CC tilt specified by a rotation matrix
   ! 30 Nov '04   add len1, len2, len3 so each edge can have a different length--D. Marx
   TYPE, PUBLIC :: CornerCube
      PRIVATE
      REAL(pr) :: size                               ! Length of an edges or scaling factor
      CHARACTER( LEN = 8 ) :: faceShape = "square"   ! ("square" or "triangle" or "sim")
      REAL(pr), DIMENSION(3) :: edgeLength = (/ 1.0_pr, 1.0_pr, 1.0_pr /) ! length of edge = size*edgeLength
      REAL(pr) :: xc = 0.0_pr                        ! Vertex location
      REAL(pr) :: yc = 0.0_pr
      REAL(pr), DIMENSION(3) :: gapWidth = (/ 0.0_pr, 0.0_pr, 0.0_pr /)
      REAL(pr) :: spin   = 0.0_pr         ! Default clocking angle
      REAL(pr), DIMENSION(3,3) :: rotate = RESHAPE( (/ 1.0_pr, 0.0_pr, 0.0_pr,    &
                                                       0.0_pr, 1.0_pr, 0.0_pr,    &
                                                       0.0_pr, 0.0_pr, 1.0_pr /), &
                                                       (/3,3/), ORDER=(/2,1/)      )
      REAL(pr) :: d1 = 0.0_pr             ! Dihedral errors (radians)
      REAL(pr) :: d2 = 0.0_pr
      REAL(pr) :: d3 = 0.0_pr
      TYPE(polygon)                     :: clearAperture
      TYPE(polygon), DIMENSION(6)       :: segPoly ! Denotes each of 6 segments
      REAL(pr),      DIMENSION(6)       :: xTiltError ! Due to dihedral error
      REAL(pr),      DIMENSION(6)       :: yTiltError
      REAL(pr), DIMENSION(3,3)          :: eM0           ! Default orientation
      REAL(pr), DIMENSION(:,:), POINTER :: fM0 => NULL() ! vertices describing face shapes
      REAL(pr), DIMENSION(3,3)          :: gap0
      REAL(pr), DIMENSION(3,3)          :: eM            ! Actual orientation
      REAL(pr), DIMENSION(:,:), POINTER :: fM  => NULL()
      REAL(pr), DIMENSION(3,3)          :: gap
   END TYPE CornerCube
   !...........................................................................

CONTAINS

  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE NewCornerCube( cc, size, shape, xc, yc,  &
                                 spin, rotMatrix,      &
                                 d1,   d2,    d3,      &
                                 gap1, gap2,  gap3,    &
                                 len1, len2,  len3      )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Initializes parameter values for a corner cube

   !     23 Oct  '00
   !     17 Apr  '01    Fortran 95
   !     18 June '01    "Rectangular" tilts: tiltX and tiltY
   !     16 Nov  '01    Back to alt-az...for Superposition Simulations
   !     15 Apr  '03    Each of 3 gaps may have separate values
   !      1 May  '03    Rotation matrix defines cc tilt orientation
   !     30 Nov  '04    add len1, len2, len3 so each edge can have a different length--D. Marx
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(CornerCube),         INTENT(IN OUT)       :: cc
   REAL(pr),                 INTENT(IN)           :: size
   CHARACTER( LEN = * ),     INTENT(IN), OPTIONAL :: shape
   REAL(pr),                 INTENT(IN), OPTIONAL :: xc, yc
   REAL(pr),                 INTENT(IN), OPTIONAL :: spin
   REAL(pr), DIMENSION(3,3), INTENT(IN), OPTIONAL :: rotMatrix
   REAL(pr),                 INTENT(IN), OPTIONAL :: d1,   d2,   d3
   REAL(pr),                 INTENT(IN), OPTIONAL :: gap1, gap2, gap3
   REAL(pr),                 INTENT(IN), OPTIONAL :: len1, len2, len3

   INTRINSIC PRESENT

   ! Local variables: none
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cc%size = size

   IF ( PRESENT( shape ) ) cc%faceShape = shape ! Shape of reflecting faces
   IF ( PRESENT( xc ) )    cc%xc = xc
   IF ( PRESENT( yc ) )    cc%yc = yc

   IF ( PRESENT( gap1 ) )    cc%gapWidth(1) = gap1
   IF ( PRESENT( gap2 ) )    cc%gapWidth(2) = gap2
   IF ( PRESENT( gap3 ) )    cc%gapWidth(3) = gap3

   IF ( PRESENT( len1 ) )  cc%edgeLength(1) = len1
   IF ( PRESENT( len2 ) )  cc%edgeLength(2) = len2
   IF ( PRESENT( len3 ) )  cc%edgeLength(3) = len3
   
   CALL Calc_eM_fM( cc )

   IF ( PRESENT( spin ) )      cc%spin   = spin    ! Orientation of corner cube
   IF ( PRESENT( rotMatrix ) ) cc%rotate = rotMatrix

   CALL Calc_Clear_Ap( cc )

   IF ( PRESENT( d1 ) )     cc%d1 = d1          ! Dihedral angle errors
   IF ( PRESENT( d2 ) )     cc%d2 = d2
   IF ( PRESENT( d3 ) )     cc%d3 = d3

   CALL Calc_Tilt_Errors( cc )                  ! Wavefront tilts

   END SUBROUTINE NewCornerCube



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Calc_eM_fM( cc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates matrices describing default orientation of corner cube. The
   ! (x,y,z) coordinate system is oriented such that the z-axis is parallel
   ! to the cube's optical axis

   !     24 Oct  '00
   !     15 Apr  '03    Gap vectors are defined
   !     01 Dec  '04    Include length of edges
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines

   IMPLICIT NONE

   TYPE(CornerCube), INTENT(IN OUT) :: cc

   INTRINSIC ASSOCIATED, ACOS, MATMUL

   ! Local variables
   REAL(pr), DIMENSION(3)   :: e1, e2, e3, c1, c2, c3
   REAL(pr), DIMENSION(3,3) :: R_x
   REAL(pr)                 :: theta, el1, el2, el3
   INTEGER                  :: NP = 10 !! number of points to describe the quarter-circle in "sim"
   INTEGER                  :: ii
   REAL(pr)                 :: edge1theta, edge1phi, psi
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   ! Unit vectors directed along edges of cube (vertex is at x = y = z = 0)
   e1 = (/  1.0_pr , 1.0_pr , 0.0_pr /)/SQRT2
   e2 = (/ -1.0_pr , 1.0_pr , 0.0_pr /)/SQRT2
   e3 = (/  0.0_pr , 0.0_pr , 1.0_pr /)

   ! length of each edge
   el1 = cc%size * cc%edgeLength(1)
   el2 = cc%size * cc%edgeLength(2)
   el3 = cc%size * cc%edgeLength(3)

   ! 
   cc%eM0(:,1) = e1
   cc%eM0(:,2) = e2
   cc%eM0(:,3) = e3

   ! Gap vectors, specifying orientation of narrow dimension of gaps
   cc%gap0(:,1) = cc%gapWidth(1)*( e2 - e3 )/SQRT2
   cc%gap0(:,2) = cc%gapWidth(2)*( e3 - e1 )/SQRT2
   cc%gap0(:,3) = cc%gapWidth(3)*( e1 - e2 )/SQRT2

   ! The eMatrix defines the orientation of the edges between faces.
   ! The fMatrix defines the shape of the effective area of the corner cube.

   SELECT CASE ( cc%faceShape )
      CASE ( "square" )
         ! Vectors directed along diagonals of square mirror faces
         c1 = el1*e1 + el2*e2
         c2 = el2*e2 + el3*e3
         c3 = el3*e3 + el1*e1
         IF ( .NOT.ASSOCIATED( cc%fM0 ) ) ALLOCATE( cc%fM0(3,6) )
         cc%fM0(:,1) = el1*e1
         cc%fM0(:,2) = c1
         cc%fM0(:,3) = el2*e2
         cc%fM0(:,4) = c2
         cc%fM0(:,5) = el3*e3
         cc%fM0(:,6) = c3
      CASE ( "triangle" )
         IF ( .NOT.ASSOCIATED( cc%fM0 ) ) ALLOCATE( cc%fM0(3,3) )
         cc%fM0(:,1) = el1*e1
         cc%fM0(:,2) = el2*e2
         cc%fM0(:,3) = el3*e3
      CASE ( "sim" )
         !! the base plate (surface normal #1) of the corner cube used in the SIM flight system
         !! is circular. The wedge mirrors form rectangles, except that the top edge has a sloped
         !! region. The height of the wedges (edge #1) is different from the length of the edges
         !! formed with the base plate. For now we will assume that all the edges have the same
         !! length (cc%size), the two faces formed by the wedges are square, and the face formed
         !! by the baseplate is a quarter circle.
         !! e1 is the edge vector normal to the baseplate. Then with a right-handed system, the 
         !! counter-clockwise list of polygon corners is e1, e1+e2, e2, points along the rounded
         !! baseplate, e3, e3+e1
         
         !! check that edge2 and edge3 have the same length, which is equal to the radius of the 
         !! baseplate of the corner cube
         IF ( abs( (el2 - el3)/el2 ) > 1.0_pr ) CALL ErrorExit( "corner cube baseplate is not circular" )

         IF ( ASSOCIATED( cc%fM0 ) ) DEALLOCATE ( cc%fM0 )
         ALLOCATE( cc%fM0(3,4+NP) )

         !! start at e1 and list the polygon points
         cc%fM0(:,1) = el1*e1
         cc%fM0(:,2) = el1*e1 + el2*e2
         !!cc%fM0(:,3) = e2  e2 will be the first point of the quarter-circle described next

         !! find the points along the quarter circle with radius = |e2|, rotated about e1
         !! the following assumes that all vectors are normalized
         edge1theta = acos(e1(3))
         edge1phi   = atan2(e1(2),e1(1))
         DO ii = 1, NP
            psi = (ii-1)*pi/2/NP

            cc%fM0(:,2+ii) = el2* &
                MATMUL(Z_Rotation(edge1phi), &
                 MATMUL(Y_Rotation(edge1theta), &
                  MATMUL(Z_Rotation(psi), &
                   MATMUL(Y_Rotation(-edge1theta), &
                    MATMUL(Z_Rotation(-edge1phi), e2 ) ))))
         ENDDO

         cc%fM0(:,2+NP+1) = el3*e3
         cc%fM0(:,2+NP+2) = el3*e3 + el1*e1

        !! OPEN(11,file="C:\tmp_vectors.txt")
        !! do ii = 1,4+NP
        !!     write(11,*) cc%fM0(:,ii)
        !! enddo
        !! Close(11)

      CASE DEFAULT
         CALL ErrorExit( "Calc_eM_fM: incorrect face shape for corner cube" )
   END SELECT

   ! Rotate so corner cube optical axis is along z (default orientation)
   theta = ACOS( 1.0_pr/SQRT3 )
   R_x = X_Rotation( theta )
   cc%eM0  = MATMUL( R_x, cc%eM0 )
   cc%fM0  = MATMUL( R_x, cc%fM0 )
   cc%gap0 = MATMUL( R_x, cc%gap0 )

   END SUBROUTINE Calc_eM_fM



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Calc_Clear_Ap( cc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Calculates a polygon describing the clear aperture of the input corner cube
   !  For corner cube geometry, see:
   !  H.D. Eckhardt, "Simple Model of Corner Reflector Phenomena,"
   !                  Appl.Opt. vol.10, 1559-1566 (1971)

   !     23 Oct  '00
   !     28 May  '01    Corrected allocation for cc%fM array
   !     18 June '01    Orientation in terms of tiltX, tiltY
   !     19 July '01    Clear polygons frontFace and backFace after use
   !     16 Nov  '01    Orientation again defined by tilt, azTilt
   !     15 Apr  '03    Gap vectors also oriented
   !      1 May  '03    Orientation defined by rotation matrix
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   TYPE(CornerCube), INTENT(IN OUT) :: cc

   INTRINSIC ASSOCIATED, SIZE, MATMUL, TINY, TRANSPOSE

   ! Local variables
   TYPE(polygon)                         :: frontFace, backFace
   REAL(pr), DIMENSION(:,:), ALLOCATABLE :: face
   REAL(pr), DIMENSION(3,3)              :: R_x, R_y, R_z, Rot
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( .NOT.ASSOCIATED( cc%fM ) ) ALLOCATE( cc%fM(3,SIZE(cc%fM0,2)) )

   ! Apply rotations to orient cube corner's optical axis
   IF ( abs( cc%spin ) > TINY(0.0_pr) )THEN
      R_z = Z_Rotation( cc%spin )                  ! Rotate about optical axis
      cc%eM  = MATMUL( R_z, cc%eM0 )
      cc%fM  = MATMUL( R_z, cc%fM0 )
   	cc%gap = MATMUL( R_z, cc%gap0 )
   ELSE
      cc%eM(:,:)  = cc%eM0(:,:)
      cc%fM(:,:)  = cc%fM0(:,:)
   	cc%gap(:,:) = cc%gap0(:,:)
   END IF
	cc%eM  = MATMUL( cc%rotate, cc%eM )
	cc%fM  = MATMUL( cc%rotate, cc%fM )
	cc%gap = MATMUL( cc%rotate, cc%gap )

   ! Polygon describing the front aperture of the cube corner, projected onto
   ! the x-y plane (containing the cube's vertex). The back face is the reflection
   ! of the front face, projected onto the x-y plane.
   ALLOCATE( face(2,SIZE(cc%fM,2)) )
   face = cc%fM(1:2,:)
   frontFace =  face
   backFace  = -face

   ! Move corner cube's vertex to its specified location
   CALL ShiftPolygon( frontFace, cc%xc, cc%yc )
   CALL ShiftPolygon( backFace,  cc%xc, cc%yc )
   cc%clearAperture = PolygonOverlap( frontFace, backFace )

   DEALLOCATE( face )
   CALL ClearPolygon( frontFace )
   CALL ClearPolygon( backFace )

   END SUBROUTINE Calc_Clear_Ap



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE Calc_Tilt_Errors( cc )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! For each of the six segments, calculates a polygon that describes its
   ! boundary, and wavefront tilt errors due to dihedral angle errors.
   !  For corner cube geometry, see:
   !  H.D. Eckhardt, "Simple Model of Corner Reflector Phenomena,"
   !                  Appl.Opt. vol.10, 1559-1566 (1971)


   !     24 Oct  '00
   !     27 May  '01    Calls MakeArray to initialize cc%segPoly()
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Array_Routines

   IMPLICIT NONE

   TYPE(CornerCube), INTENT(IN OUT) :: cc

   INTRINSIC ABS, TINY, TRANSPOSE, MATMUL, DOT_PRODUCT

   ! Local variables
   REAL(pr), DIMENSION(2,7) :: segVertex
   REAL(pr), DIMENSION(3,3) :: Rx, Ry, Rz, RxT, RyT, RzT, M
   REAL(pr), DIMENSION(3)   :: xUnitVector, yUnitVector, zUnitVector
   REAL(pr), DIMENSION(3)   :: retroDir
   INTEGER                  :: k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IF ( ( ABS( cc%d1 ) <= TINY(0.0_pr) ) .AND.  &
        ( ABS( cc%d2 ) <= TINY(0.0_pr) ) .AND.  &
        ( ABS( cc%d3 ) <= TINY(0.0_pr) ) ) THEN
      DO k = 1,6
         cc%xTiltError(k) = 0.0_pr
         cc%yTiltError(k) = 0.0_pr
      END DO
      RETURN
   END IF

   ! Calculate tilt errors in each segment, due to dihedral angle errors.
   ! The following vectors describe the boundaries between segments:

   segVertex(:,1) =  (2.0_pr*cc%size*cc%edgeLength(1))*cc%eM(1:2,1)
   segVertex(:,2) = -(2.0_pr*cc%size*cc%edgeLength(3))*cc%eM(1:2,3)
   segVertex(:,3) =  (2.0_pr*cc%size*cc%edgeLength(2))*cc%eM(1:2,2)
   segVertex(:,4) = -(2.0_pr*cc%size*cc%edgeLength(1))*cc%eM(1:2,1)
   segVertex(:,5) =  (2.0_pr*cc%size*cc%edgeLength(3))*cc%eM(1:2,3)
   segVertex(:,6) = -(2.0_pr*cc%size*cc%edgeLength(2))*cc%eM(1:2,2)
   segVertex(:,7) = segVertex(:,1)
   DO k = 1,7
      segVertex(1,k) = segVertex(1,k) + cc%xc
      segVertex(2,k) = segVertex(2,k) + cc%yc
   END DO

   Rx = X_Rotation( -2.0_pr * cc%d1 )  ! Clockwise rotation for positive error
   Ry = Y_Rotation( -2.0_pr * cc%d2 )
   Rz = Z_Rotation( -2.0_pr * cc%d3 )
   RxT = TRANSPOSE( Rx )
   RyT = TRANSPOSE( Ry )
   RzT = TRANSPOSE( Rz )

   ! Unit vectors for (x,y,z) expressed in corner cube coordinates (the incident
   ! wavefront amplitude will be specified in the x-y plane, with propagation
   ! along the minus z-axis).
   xUnitVector = cc%eM(1,:)
   yUnitVector = cc%eM(2,:)
   zUnitVector = cc%eM(3,:)

   DO k = 1,6
      ! Create a polygon (triangle) delimiting each of the 6 segments
      cc%segPoly(k) = MakeArray( 2, 3,                                         &
      							 (/ cc%xc, segVertex(1,k), segVertex(1,k+1),   &
                                    cc%yc, segVertex(2,k), segVertex(2,k+1) /) )
      ! Calculate tilt for this segment -- see Eckhardt (1971)
      SELECT CASE ( k )
         case (1)
            M = MATMUL( Rx, MATMUL( Rz, RyT ) )
         case (2)
            M = MATMUL( Rx, MATMUL( RyT, RzT ) )
         case (3)
            M = MATMUL( Ry, MATMUL( Rx, RzT ) )
         case (4)
            M = MATMUL( Ry, MATMUL( RzT, RxT ) )
         case (5)
            M = MATMUL( Rz, MATMUL( Ry, RxT ) )
         case (6)
            M = MATMUL( Rz, MATMUL( RxT, RyT ) )
         CASE DEFAULT
            CALL ErrorExit( "Calc_Tilt_Errors: " // &
                        "incorrect number of corner cube segments" )
      END SELECT
      ! Cube's optical axis in wavefront (x,y,z) coordinates
      retroDir = MATMUL( M, zUnitVector )

      cc%xTiltError(k) = DOT_PRODUCT( retroDir, xUnitVector )
      cc%yTiltError(k) = DOT_PRODUCT( retroDir, yUnitVector )
   END DO

   END SUBROUTINE Calc_Tilt_Errors



  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SetClockAngle( cc, spinAngle )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Rotates corner cube about its optical axis; this routine and SetTilt
   ! replace former routine, SetOrientation()

   !      1 May '03
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(CornerCube), INTENT(IN OUT) :: cc
   REAL(pr),         INTENT(IN)     :: spinAngle

   ! Local variables: none
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cc%spin   = spinAngle

   CALL Calc_Clear_Ap( cc )
   CALL Calc_Tilt_Errors( cc )

   END SUBROUTINE SetClockAngle


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SetTilt( cc, rotMatrix )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Tilts corner cube; this routine and SetClockAngle replace former routine,
   ! SetOrientation()

   !      1 May '03
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(CornerCube),     INTENT(IN OUT) :: cc
   REAL(pr), DIMENSION(3,3), INTENT(IN) :: rotMatrix

   ! Local variables: none
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cc%rotate = rotMatrix

   CALL Calc_Clear_Ap( cc )
   CALL Calc_Tilt_Errors( cc )

   END SUBROUTINE SetTilt




  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE SetDihedralErrors( cc, d1, d2, d3 )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Assigns values to corner cube dihedral errors and calculates corresponding
   ! wavefront tilt errors

   !     24 Oct  '00
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(CornerCube), INTENT(IN OUT) :: cc
   REAL(pr),         INTENT(IN)     :: d1, d2, d3

   ! Local variables: none
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cc%d1 = d1
   cc%d2 = d2
   cc%d3 = d3

   CALL Calc_Tilt_Errors( cc )

   END SUBROUTINE SetDihedralErrors


  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE ShiftCornerCubeVertex( cc, xShift, yShift )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   ! Shifts the location of the vertex
   !     27 Feb '02
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   IMPLICIT NONE

   TYPE(CornerCube), INTENT(IN OUT) :: cc
   REAL(pr),         INTENT(IN)     :: xShift, yShift

   ! Local variables:
   INTEGER :: k
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   cc%xc = cc%xc + xShift
   cc%yc = cc%yc + yShift

   CALL ShiftPolygon( cc%clearAperture, xShift, yShift )

   IF ( ( ABS( cc%d1 ) > EPSILON(0.0_pr) ) .OR.  &
        ( ABS( cc%d2 ) > EPSILON(0.0_pr) ) .OR.  &
        ( ABS( cc%d3 ) > EPSILON(0.0_pr) ) )  THEN
      DO k = 1,6
   	   CALL ShiftPolygon( cc%segPoly(k), xShift, yShift )
      END DO
   END IF

   END SUBROUTINE ShiftCornerCubeVertex




  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
   SUBROUTINE CornerCubeReflect( cc, amp, dx, x, dy, y, wavelength )
  != = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =

   !  Inverts input wave amplitude through corner cube vertex, clips
   !  by corner cube clear aperture, and applies OPD due to dihedral errors.

   !  	cc    		Corner cube from which beam is reflected
   !  	amp   		Beam's complex wave amplitude
   !     dx, dy   	Sample spacing
   !  	x,  y			Coordinates of beam amplitude samples
   !     wavelength  Beam wavelength

   !     26 Oct  '00
   !     27 May  '01    Corrected testing of d1, d2, d3
   !     30 Nov  '01    Use DO loops for adding OPD and segOPD arrays
   !      3 Dec  '01    OPD and gaps are "out" only
   !      4 Dec  '01    Rewritten
   !     15 Apr  '03    Each gap may have separate width value; includes
   !                    variation of projected gap width with cc orientation
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   USE Constants
   USE Array_Routines
   USE Clipping_Routines

   IMPLICIT NONE

   TYPE(CornerCube),             INTENT(IN)    :: cc
   COMPLEX(prc), DIMENSION(:,:), INTENT(INOUT) :: amp
   REAL(pr),                     INTENT(IN)    :: dx, dy
   REAL(pr),     DIMENSION(:),   INTENT(IN)    ::  x,  y
   REAL(pr),                     INTENT(IN)    :: wavelength


   INTRINSIC EPSILON, ABS, ATAN2, SQRT, SIZE, CMPLX

   ! Local variables
   REAL(pr), DIMENSION(:,:), ALLOCATABLE :: segOPD, OPD
   REAL(pr)                             :: const
   REAL(pr)                             :: edgeAngle, gapLength, projectedWidth
   TYPE(polygon)                        :: gap
   INTEGER                              :: Nx, Ny
   INTEGER                              :: j,  k, seg, edge
   !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   open(11,file="C:\edge_vectors.dat")
   !! write out three column vectors, each vector is an edge direction
   do j = 1,3
       write(11,*) cc%eM(j,1), cc%eM(j,2), cc%eM(j,3)
   enddo

   CALL ReflectThroughOrigin( amp, cc%xc/dx, cc%yc/dy )
   CALL PolygonClip( amp, x, y, cc%clearAperture )

   IF ( ( ABS( cc%d1 ) > EPSILON(0.0_pr) ) .OR.  &
        ( ABS( cc%d2 ) > EPSILON(0.0_pr) ) .OR.  &
        ( ABS( cc%d3 ) > EPSILON(0.0_pr) ) )  THEN

      Nx = SIZE( x )
      Ny = SIZE( y )

      ALLOCATE( segOPD(Ny,Nx), OPD(Ny,Nx) )

      ! Calculate OPD and apply to amplitude; within each segment,
      ! the OPD is linear in x and y (pure tilt).

      DO k = 1,Nx
         DO j = 1,Ny
      		OPD(j,k) = 0.0_pr
         END DO
      END DO

      DO seg = 1,6

         DO k = 1,Nx
         	DO j = 1,Ny
               segOPD(j,k) = cc%xTiltError(seg) * x(k) +  &
                             cc%yTiltError(seg) * y(j)
            END DO
         END DO
         CALL PolygonClip( segOPD, x, y, cc%segPoly(seg) )
         DO k = 1,Nx
         	DO j = 1,Ny
         		OPD(j,k) = OPD(j,k) + segOPD(j,k)
         	END DO
         END DO

      END DO

      const = TWOPI/wavelength
      DO k = 1,Nx
      	DO j = 1, Ny
         	amp(j,k) = amp(j,k) * EXP( CMPLX( 0.0_pr, const*OPD(j,k), KIND=prc ) )
      	END DO
      END DO

      ! clear variables OPD and segOPD
      DEALLOCATE ( segOPD, OPD )

   END IF ! if dihedral

   ! Generate obscurations due to gaps, if any

   DO edge = 1,3

        gapLength = 2.0_pr * cc%size * cc%edgeLength(edge)
        
		IF ( cc%gapWidth(edge) > EPSILON(0.0_pr) ) THEN
			! Find orientation of corner cube edge, projected onto x-y plane
            edgeAngle = ATAN2( cc%eM(2,edge), cc%eM(1,edge) )
            projectedWidth = SQRT( cc%gap(1,edge)**2 + cc%gap(2,edge)**2 )
            gap = MakeRectangle( gapLength, projectedWidth,  &
                              edgeAngle, cc%xc, cc%yc      )
            CALL PolygonObsc( amp, x, y, gap )
            CALL ClearPolygon( gap )
        END IF
   
   END DO

   END SUBROUTINE CornerCubeReflect

END MODULE Corner_Cubes
