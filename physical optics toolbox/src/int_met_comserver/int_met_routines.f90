MODULE IntMetRoutines

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Propagates to lens, Applies lens and propagates one focal length behind lens


   USE Kinds
   USE SI_Units
   USE Constants
   USE Wavefronts
   USE Optics_Routines
   USE Utility_Routines
   USE Error_Exit


   IMPLICIT NONE

   PRIVATE


    !! AVAILABLE SUBROUTINES:
    PUBLIC :: FocusLens
    PUBLIC :: ApplyMask
    PUBLIC :: ApplyMaskRotate
    PUBLIC :: ApplyMaskGeneral
    PUBLIC :: ApplyMaskX
    PUBLIC :: ApplyMaskXRounded
    PUBLIC :: ApplyMaskPoly

   !! unit #'s for outputs
   INTEGER, PARAMETER :: iUNF        = 1   ! ionumber for writing UNF outputs


CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE FocusLens( beam, focuslens_f, focuslens_D, dxout, dyout )

IMPLICIT NONE

   type(Wavefront), INTENT(INOUT) :: beam  ! the beam
   REAL(pr)       , INTENT(IN)    :: focuslens_f ! focal length
   REAL(pr)       , INTENT(IN)    :: focuslens_D ! lens diameter
   REAL(pr)       , INTENT(IN)    :: dxout, dyout ! desired grid spacing at focal plane

   TYPE(simpleLens)    :: hetfocuslens              ! lens object
   REAL(pr)            :: dx, dy
   INTEGER             :: Nx, Ny

   CALL GetWavefrontSampling( beam, Nx, Ny, dx, dy )


    !! propagate to the next element (focusing lens)
    call Propagate( beam, focuslens_f, dx, dy, applyCurv=.TRUE. )

   !! initialize focus lens and apply focusing lens as a subroutine
   call SetLensProperties( hetfocuslens, focuslens_f, focuslens_D )

   call PropagateToFocalPlane( hetfocuslens, beam, dxout, dyout, applyCurv=.TRUE. )


END SUBROUTINE FocusLens

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine for applying mask
!! creates two rectangles, either east-west or north-south, depending on direction parameter
!! len_r = length of each rectangle in y-direction for 'ns', or x-direction for 'ew'
!! len_t = length of each rectangel in x-direction for 'ns', or y-direction for 'ns'
!! offset = distance from origin to center of each rectangle
!! direction = 'ns' 'NS' or 'ew' 'EW'
!! rotangle (optional) = angle of rotation of rectangle about its center, i.e. rotation is
!!    before translation
!! tiltAngle, tiltOrientation tilt the whole mask with respect to the propagation axis (Z). The tilt
!!    is applied last. The tiltOrientation is the angle between the x-axis and the axis of rotation.
SUBROUTINE ApplyMask(beam,len_r,len_t,offset,direction,rotangle,x_align,y_align,tiltAngle,tiltOrientation)

IMPLICIT NONE

    type(Wavefront),  INTENT(INOUT)        :: beam
    REAL(pr),         INTENT(IN)           :: len_r,len_t,offset
    CHARACTER(LEN=*), INTENT(IN)           :: direction ! "ns" or "ew"
    REAL(pr),         INTENT(IN), OPTIONAL :: rotangle ! rotates mask around origin
    REAL(pr),         INTENT(IN), OPTIONAL :: x_align, y_align  ! misalignment offset in x and y
    REAL(pr),         INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

 
    !! local variables
    type(Wavefront) :: tmpbeam
    REAL(pr)        :: len_x, len_y, off_x, off_y, phi, xc, yc, tiltA, tiltO

    !! determine offset is north-south or east-west
    SELECT CASE ( direction )
        CASE( "ns", "NS", "north-south", "n2" )
            print*, "applying mask for north-south"
            off_x = 0.0_pr
            off_y = offset
            len_x = len_t
            len_y = len_r
        CASE( "ew", "EW", "east-west", "n1" )
            print*, "applying mask for east-west"
            off_x = offset
            off_y = 0.0_pr
            len_x = len_r
            len_y = len_t
        CASE DEFAULT
            call ErrorExit( "direction must be n1 or n2" )
    END SELECT

    !! check rotation angle
    IF ( PRESENT( rotangle ) ) THEN
        phi = rotangle
    ELSE
        phi = 0.0_pr
    ENDIF

    !! check for misalignment parameters
    xc = 0.0_pr
    yc = 0.0_pr
    IF (PRESENT( x_align ) ) THEN
        xc = x_align
    ENDIF
    IF (PRESENT( y_align ) ) THEN
        yc = y_align
    ENDIF

    !! check for mask tilt
    tiltA = 0.0_pr
    tiltO = 0.0_pr
    IF (PRESENT(tiltAngle) .AND. PRESENT(tiltOrientation)) THEN
        tiltA = tiltAngle
        tiltO = tiltOrientation
    ENDIF
        
    tmpbeam = beam
    call ClipRect(tmpbeam,len_x,len_y, off_x+xc, off_y+yc, angle=phi, tiltAngle=tiltA, tiltOrientation=tiltO)
    call ClipRect(beam,   len_x,len_y,-off_x+xc,-off_y+yc, angle=phi, tiltAngle=tiltA, tiltOrientation=tiltO)

    call WavefrontCumSum(beam,tmpbeam)

    call ClearWavefront(tmpbeam)

END SUBROUTINE ApplyMask

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine for applying a rotated and tilted mask
!! creates two rectangles, either east-west or north-south, depending on direction parameter,
!! and rotates the entire mask around the origin
!! len_r = length of each rectangle in y-direction for 'ns', or x-direction for 'ew'
!! len_t = length of each rectangel in x-direction for 'ns', or y-direction for 'ns'
!! offset = distance from origin to center of each rectangle
!! rotangle = angle of rotation of the mask about the origin, i.e rotation is after the 
!!    rectangles are translated
!! direction = 'ns' 'NS' 'n2' or 'ew' 'EW' 'n1'
!! x_align, y_align are translation misalignment of the whole mask. These decenters are applied
!!    after the mask is assembled and rotated, but before tilted.
!! tiltAngle, tiltOrientation is the tilt of the mask plane w.r.t. the axis of propagation.
!!    tiltOrientation is the angle between the axis of tilt (normal to z) and the x-axis
SUBROUTINE ApplyMaskRotate(beam,len_r,len_t,offset,rotangle,direction,x_align,y_align, &
                           tiltAngle, tiltOrientation)

use Polygons

IMPLICIT NONE

    type(Wavefront),  INTENT(INOUT)        :: beam
    REAL(pr),         INTENT(IN)           :: len_r,len_t,offset
    REAL(pr),         INTENT(IN)           :: rotangle ! rotates mask around origin
    CHARACTER(LEN=*), INTENT(IN)           :: direction ! "ns" or "ew"
    REAL(pr),         INTENT(IN), OPTIONAL :: x_align, y_align  ! misalignment offset in x and y
    REAL(pr),         INTENT(IN), OPTIONAL :: tiltAngle, tiltOrientation

 
    !! local variables
    type(Wavefront)          :: tmpbeam
    REAL(pr), DIMENSION(2,4) :: vu, vd
    type(polygon)            :: rectu, rectd
    REAL(pr)                 :: xc, yc
    REAL(pr)                 :: tiltA, tiltO

    !! determine offset is north-south or east-west
    !! determine the corners of the rectangles
    !! rotate
    !! use polygon clip
    SELECT CASE ( direction )
        CASE( "ns", "NS", "north-south", "n2" )
            print*, "applying mask for north-south"
            vu(1,:) = 0.5_pr*(/-len_t,  len_t,  len_t, -len_t/)
            vu(2,:) = 0.5_pr*(/-len_r, -len_r,  len_r,  len_r/) + offset
            vd(1,:) = vu(1,:)
            vd(2,:) = 0.5_pr*(/-len_r, -len_r,  len_r,  len_r/) - offset
        CASE( "ew", "EW", "east-west", "n1" )
            print*, "applying mask for east-west"
            vu(1,:) = 0.5_pr*(/-len_r,  len_r,  len_r, -len_r/) + offset
            vu(2,:) = 0.5_pr*(/-len_t, -len_t,  len_t,  len_t/)
            vd(1,:) = 0.5_pr*(/-len_r,  len_r,  len_r, -len_r/) - offset
            vd(2,:) = vu(2,:)
        CASE DEFAULT
            call ErrorExit( "direction must be ns or ew" )
    END SELECT
    rectu = vu
    rectd = vd

    ! check for presence of x_align, y_align
    xc = 0.0_pr
    yc = 0.0_pr
    IF (PRESENT( x_align ) ) THEN
        xc = x_align
    ENDIF
    IF (PRESENT( y_align ) ) THEN
        yc = y_align
    ENDIF

    ! check for presence of tilt
    tiltA = 0.0_pr
    tiltO = 0.0_pr
    IF (PRESENT( tiltAngle ) ) THEN
        tiltA = tiltAngle
    ENDIF
    IF (PRESENT( tiltOrientation ) ) THEN
        tiltO = tiltOrientation
    ENDIF


    tmpbeam = beam
    call ClipPoly(tmpbeam, rectu, xc, yc, rotangle, &
                  tiltAngle=tiltA, tiltOrientation=tiltO)
    call ClipPoly(beam,    rectd, xc, yc, rotangle, &
                  tiltAngle=tiltA, tiltOrientation=tiltO)

    call WavefrontCumSum(beam,tmpbeam)

    call ClearWavefront(tmpbeam)

END SUBROUTINE ApplyMaskRotate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine for applying a more general mask consisting of two rectangles symmetric about the origin
SUBROUTINE ApplyMaskGeneral(beam, len_x, len_y, off_x, off_y)

IMPLICIT NONE

    type(Wavefront), INTENT(INOUT) :: beam
    REAL(pr), INTENT(IN)           :: len_x, len_y, off_x, off_y

    type(Wavefront) :: tmpbeam

    tmpbeam = beam
    call ClipRect(tmpbeam,len_x,len_y, off_x, off_y)
    call ClipRect(beam,   len_x,len_y,-off_x,-off_y)

    call WavefrontCumSum(beam,tmpbeam)

    call ClearWavefront(tmpbeam)

END SUBROUTINE ApplyMaskGeneral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Subroutine for applying a mask described by polygon vertices
!! Applies the ClipPolygon twice, once for the given polygon, and once for the
!! polygon reflected about the origin (x => -x, y => -y)
!! vertices is a 2 row x N column list of x,y pairs of vertices of the polygon
!! the vertices must be in counter-clockwise order
!! xc, yc (optional) are x and y offsets (misalignments) of the mask
SUBROUTINE ApplyMaskPoly(beam, vertices, xc, yc)

use Polygons

IMPLICIT NONE

    type(Wavefront),            INTENT(INOUT)  :: beam
    REAL(pr), DIMENSION(1:,1:), INTENT(IN)     :: vertices
    REAL(pr),                   INTENT(IN)     :: xc, yc
 
    !! local variables
    type(Wavefront)          :: tmpbeam
    type(polygon)            :: tmppoly

    tmpbeam = beam

    tmppoly = vertices
    call ClipPoly(tmpbeam, tmppoly, xc, yc)

    tmppoly = -vertices
    call ClipPoly(beam,    tmppoly, xc, yc)

    call WavefrontCumSum(beam,tmpbeam)

    call ClearWavefront(tmpbeam)

END SUBROUTINE ApplyMaskPoly

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Mask consisting of four squares set on an 'X' pattern
!! 
!! radius = distance from origin to center of each of the four squares
!! length = side of each square
SUBROUTINE ApplyMaskX(beam, radius, length)

IMPLICIT NONE

    type(Wavefront), INTENT(INOUT) :: beam
    REAL(pr), INTENT(IN)           :: radius, length

    type(Wavefront) :: tmpa, tmpb
    REAL(pr)        :: xc, yc

    ! calculate xc, yc offsets from radius
    xc = radius/SQRT(2.0_pr)
    yc = radius/SQRT(2.0_pr)

    tmpa = beam
    tmpb = beam
    call ClipRect(beam, length, length, xc, yc)
    call ClipRect(tmpa, length, length,-xc, yc)
    call WavefrontCumSum(beam,tmpa)

    tmpa = tmpb
    call ClipRect(tmpa, length, length, xc,-yc)
    call WavefrontCumSum(beam,tmpa)

    call ClipRect(tmpb, length, length,-xc,-yc)
    call WavefrontCumSum(beam,tmpb)

    call ClearWavefront(tmpa)
    call ClearWavefront(tmpb)

END SUBROUTINE ApplyMaskX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Mask consisting of four squares set on an 'X' patter, where the corner closest to the
!! of each square is rounded using a specified radius
!!
!! radius = distance from origin to center of each of the four squares
!! length = side of each square
!! corner_rad = radius of curvature of the rounded corners
SUBROUTINE ApplyMaskXRounded(beam, radius, length, corner_rad)

    USE Polygons

IMPLICIT NONE

    type(Wavefront), INTENT(INOUT) :: beam
    REAL(pr), INTENT(IN)           :: radius, length, corner_rad

    type(Wavefront)    :: tmpa
    INTEGER            :: Nx, Ny
    REAL(pr)           :: dx, dy
    REAL(pr)           :: lambda

    REAL(pr)           :: xc, yc    ! center of square hole
    REAL(pr)           :: xr, yr    ! center of circle for rounded corner
    REAL(pr)           :: xd, yd    ! center offsets for rectangles

    ! calculate xc, yc offsets from radius
    xc = radius/SQRT(2.0_pr)           ! center of the square
    yc = radius/SQRT(2.0_pr)
    xr = xc - length/2 + corner_rad    ! center of the circle
    yr = yc - length/2 + corner_rad
    xd = xc + corner_rad/2             ! center offset for rectangles
    yd = yc + corner_rad/2

    ! create a negative mask using obscurations
    ! first create a wavefront = 1
    call GetWavefrontSampling( beam, Nx, Ny, dx, dy )
    lambda = WavefrontWavelength( beam )
    call GenTopHatBeam( tmpa, Nx, Ny, dx, dy, lambda, 2*(Nx*dx + Ny*dy) )

    ! obscure the four circles
    call ClipCirc( tmpa, 2*corner_rad, xr, yr, obsc = .TRUE. )
    call ClipCirc( tmpa, 2*corner_rad,-xr, yr, obsc = .TRUE. )
    call ClipCirc( tmpa, 2*corner_rad, xr,-yr, obsc = .TRUE. )
    call ClipCirc( tmpa, 2*corner_rad,-xr,-yr, obsc = .TRUE. )

    ! obscure rectangular regions to make the rest of the square
    call ClipRect( tmpa, length-corner_rad, length, xd, yc, obsc = .TRUE. )
    call ClipRect( tmpa, length, length-corner_rad, xc, yd, obsc = .TRUE. )

    call ClipRect( tmpa, length-corner_rad, length,-xd, yc, obsc = .TRUE. )
    call ClipRect( tmpa, length, length-corner_rad,-xc, yd, obsc = .TRUE. )

    call ClipRect( tmpa, length-corner_rad, length, xd,-yc, obsc = .TRUE. )
    call ClipRect( tmpa, length, length-corner_rad, xc,-yd, obsc = .TRUE. )

    call ClipRect( tmpa, length-corner_rad, length,-xd,-yc, obsc = .TRUE. )
    call ClipRect( tmpa, length, length-corner_rad,-xc,-yd, obsc = .TRUE. )

    ! apply the negative mask to the beam
    call ApplyArbitraryMask( beam, tmpa, negative = .TRUE. )


END SUBROUTINE ApplyMaskXRounded

END MODULE IntMetRoutines
