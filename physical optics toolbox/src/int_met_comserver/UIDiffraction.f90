!
!  UIDiffraction.f90 - This file contains the implementation of the
!                    IDiffraction methods
!

    ! IDiffraction_x_vector_get
    function IDiffraction_x_vector_get( ObjectData ,&
             BeamID,&
			 x_vector) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamID
        REAL(8), intent(inout) :: x_vector
        DIMENSION x_vector(1:,1:)
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! local variables
        integer nxin, nyin, nxwv
        REAL(pr),     ALLOCATABLE, DIMENSION(:)   :: wvfront_x

        ! check for valid BeamID and
        ! get the x-coord vector from the wavefront and get its size
        IF ( BeamID >= 0 .and. BeamID <= MAX_NR_WFRONTS) THEN
            call GetWavefrontCoords(ObjectData%Wbeams(BeamID), x = wvfront_x)

        ELSE
            ! BeamID error
            allocate(wvfront_x(1))
            wvfront_x(1) = 0.0_pr
        ENDIF

        ! size of output data
        nxwv = SIZE(wvfront_x)

        ! get nx of the client supplied vector
        nxin = SIZE(x_vector,1)
        nyin = SIZE(x_vector,2)

        ! check which dimension is for the vector (Matlab always uses 2-d matrix,
        ! even for vectors)
        IF (nxin == nxwv) THEN
            x_vector(:,1) = wvfront_x
        ELSEIF (nyin == nxwv) THEN
            x_vector(1,:) = wvfront_x
        ELSE 
            ! error condition
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_y_vector_get
    function IDiffraction_y_vector_get( ObjectData ,&
             BeamID,&
			 y_vector) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamID
        REAL(8), intent(inout) :: y_vector
        DIMENSION y_vector(1:,1:)
        integer(LONG) hresult
        ! TODO:  Add implementation

        integer nxin, nyin, nywv
        REAL(pr),     ALLOCATABLE, DIMENSION(:)   :: wvfront_y

        ! check for valid BeamID and
        ! get the y-coord vector from the wavefront and get its size
        IF ( BeamID >= 0 .and. BeamID <= MAX_NR_WFRONTS) THEN
            call GetWavefrontCoords(ObjectData%Wbeams(BeamID), y = wvfront_y)

        ELSE
            ! BeamID error
            allocate(wvfront_y(1))
            wvfront_y(1) = 0.0_pr
        ENDIF

        ! size of output data
        nywv = SIZE(wvfront_y)

        ! get nx of the client supplied vector
        nxin = SIZE(y_vector,1)
        nyin = SIZE(y_vector,2)

        ! check which dimension is for the vector (Matlab always uses 2-d matrix,
        ! even for vectors)
        IF (nxin == nywv) THEN
            y_vector(:,1) = wvfront_y
        ELSEIF (nyin == nywv) THEN
            y_vector(1,:) = wvfront_y
        ELSE 
            ! error condition
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_wavefront_put
    function IDiffraction_wavefront_put( ObjectData ,&
             amp_r,&
			 amp_i,&
			 dx,&
			 dy,&
			 wavelength,&
			 curv,&
			 BeamID,&
			 status) result (hresult)
        use Diffraction_Types
        use KINDS
        use Wavefronts
        use Error_Exit

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: amp_r
        DIMENSION amp_r(1:,1:)
        REAL(8), intent(in) :: amp_i
        DIMENSION amp_i(1:,1:)
        REAL(8), intent(in) :: dx
        REAL(8), intent(in) :: dy
        REAL(8), intent(in) :: wavelength
        REAL(8), intent(in) :: curv
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
        
        ! local variables
        ! temporary complex array to hold amp_r + j*amp_i
        COMPLEX(prc), ALLOCATABLE, DIMENSION(:,:) :: tmpbeam
        INTEGER                                   :: ny, nx

        ! initialize status to no error
        status = INTMETCOMNOERROR

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE

            ! get size of wavefront Ny,Nx and check that amp_r and amp_i agree
            IF ( SIZE(amp_r,1).EQ.SIZE(amp_i,1) .AND. &
                 SIZE(amp_r,2).EQ.SIZE(amp_i,2)) THEN
    
                 ny = SIZE(amp_r,1)
                 nx = SIZE(amp_r,2)
    
                 ALLOCATE(tmpbeam(ny,nx))
                 
                 tmpbeam = amp_r + CMPLX(0.0_pr,1.0_pr,prc)*amp_i

                 call CreateWavefront(ObjectData%Wbeams(BeamID), &
                    tmpbeam, dx, dy, wavelength, curv)

                 DEALLOCATE(tmpbeam)

            ELSE
                status = INVALIDWVSIZE

            ENDIF


        ENDIF ! if valid BeamID

        hresult = S_OK
    end function

    ! IDiffraction_wavefront_get
    function IDiffraction_wavefront_get( ObjectData ,&
             BeamID,&
			 amp_r,&
			 amp_i) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamID
        REAL(8), intent(inout) :: amp_r
        DIMENSION amp_r(1:,1:)
        REAL(8), intent(inout) :: amp_i
        DIMENSION amp_i(1:,1:)
        integer(LONG) hresult
        ! TODO:  Add implementation
        ! BeamID = 0 for standard wavefront

         ! local variables
        COMPLEX(prc), ALLOCATABLE, DIMENSION(:,:) :: tmpbeam
        INTEGER                                   :: nx, ny

        ! check BeamID
        IF ( BeamID >= 0 .and. BeamID <= MAX_NR_WFRONTS ) THEN
            call GetWavefrontSampling(ObjectData%Wbeams(BeamID), Nx=nx, Ny=ny)
            ALLOCATE(tmpbeam(ny,nx))

            tmpbeam = GetWavefrontAmp(ObjectData%Wbeams(BeamID), nx, ny)
            
        ELSE
            ! BeamID error
            ALLOCATE(tmpbeam(1,1))
            tmpbeam(1,1) = 0.0_pr
        
        ENDIF

        amp_r = REAL(tmpbeam)
        amp_i = IMAG(tmpbeam)

        ! de-allocate temporary array
        DEALLOCATE(tmpbeam)

        hresult = S_OK
    end function

    ! IDiffraction_CopyWavefront
    function IDiffraction_CopyWavefront( ObjectData ,&
             BeamIDfrom,&
			 BeamIDto,&
			 status) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamIDfrom
        INTEGER(4), intent(in) :: BeamIDto
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! initialize to no error
        status = INTMETCOMNOERROR

        ! BeamID = 0 for standard wavefront
        IF ( BeamIDfrom .eq. BeamIDto ) THEN
            ! do nothing

        ELSEIF ( BeamIDfrom >= 0 .AND. BeamIDfrom <= MAX_NR_WFRONTS .AND. &
                 BeamIDto   >= 0 .AND. BeamIDto   <= MAX_NR_WFRONTS ) THEN
            ObjectData%Wbeams(BeamIDto) = ObjectData%Wbeams(BeamIDfrom)

        ELSE
            ! error condition
            status = INVALIDWFRONTID
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_CopyWavefrontAmp
    ! replace the amplitude of BeamIDPhase with the amplitude of BeamIDAmp
    ! BeamIDPhase is overwritten, BeamIDAmp is not.
    function IDiffraction_CopyWavefrontAmp( ObjectData ,&
             BeamIDAmp,&
			 BeamIDPhase,&
			 status) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamIDAmp
        INTEGER(4), intent(in) :: BeamIDPhase
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! initialize to no error
        status = INTMETCOMNOERROR

        ! if BeamID's are valid
        IF ( BeamIDAmp .eq. BeamIDPhase ) THEN
            ! do nothing

        ELSEIF ( BeamIDAmp   >= 0 .AND. BeamIDAmp   <= MAX_NR_WFRONTS .AND. &
                 BeamIDPhase >= 0 .AND. BeamIDPhase <= MAX_NR_WFRONTS ) THEN

            call CopyWavefrontAmp(ObjectData%Wbeams(BeamIDAmp), &
                                  ObjectData%Wbeams(BeamIDPhase) )

        ELSE
            ! error condition
            status = INVALIDWFRONTID
        ENDIF


        hresult = S_OK
    end function

    ! IDiffraction_SumWavefront
    function IDiffraction_SumWavefront( ObjectData ,&
             BeamID,&
			 BeamIDadd,&
			 status) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(in) :: BeamIDadd
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! initialize to no error
        status = INTMETCOMNOERROR

        ! BeamID = 0 for standard wavefront
        IF ( BeamID .eq. BeamIDadd ) THEN
            ! cannot add a wavefront to itself
            status = INVALIDWFRONTID

        ELSEIF ( BeamID    >= 0 .AND. BeamID    <= MAX_NR_WFRONTS .AND. &
                 BeamIDadd >= 0 .AND. BeamIDadd <= MAX_NR_WFRONTS ) THEN

            call WavefrontCumSum(ObjectData%Wbeams(BeamID),ObjectData%Wbeams(BeamIDadd))

        ELSE
            ! error condition
            status = INVALIDWFRONTID
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_GetWavefrontParms
    function IDiffraction_GetWavefrontParms( ObjectData ,&
             BeamID,&
			 Nx,&
			 Ny,&
			 dx,&
			 dy,&
			 curv,&
			 wavelength) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: Nx
        INTEGER(4), intent(out) :: Ny
        REAL(8), intent(out) :: dx
        REAL(8), intent(out) :: dy
        REAL(8), intent(out) :: curv
        REAL(8), intent(out) :: wavelength
        integer(LONG) hresult
        ! TODO:  Add implementation

        
        IF ( BeamID >= 0 .AND. BeamID <= MAX_NR_WFRONTS ) THEN

            call GetWavefrontSampling(ObjectData%Wbeams(BeamID), Nx, Ny, dx, dy)
            curv = WavefrontCurvature(ObjectData%Wbeams(BeamID))
            wavelength = WavefrontWavelength(ObjectData%Wbeams(BeamID))

        ELSE
            ! BeamID error
            Nx = 0
            Ny = 0
            dx = 0.0_pr
            dy = 0.0_pr
            curv = 0.0_pr
            wavelength = 0.0_pr

        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_WriteWavefrontUNF
    function IDiffraction_WriteWavefrontUNF( ObjectData ,&
             StrFilename,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        CHARACTER(*), intent(in) :: StrFilename
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE

            call WriteWavefrontUNF( ObjectData%Wbeams(BeamID), StrFilename )

            status = INTMETCOMNOERROR

        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ReadWavefrontUNF
    function IDiffraction_ReadWavefrontUNF( ObjectData ,&
             StrFilename,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use Wavefronts

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        CHARACTER(*), intent(in) :: StrFilename
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE

            call ReadWavefrontUNF( ObjectData%Wbeams(BeamID), StrFilename )

            status = INTMETCOMNOERROR

        ENDIF

                
        hresult = S_OK
    end function

    ! IDiffraction_WavefrontPower
    function IDiffraction_WavefrontPower( ObjectData ,&
             power,&
			 BeamID) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(out) :: power
        INTEGER(4), intent(in) :: BeamID
        integer(LONG) hresult
        ! TODO:  Add implementation
        
        IF ( BeamID >= 0 .and. BeamID <= MAX_NR_WFRONTS ) THEN
            power = WavefrontPower(ObjectData%Wbeams(BeamID))
        ELSE
            ! BeamID error
            power = 0.0_pr
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_CreateGausSource
    function IDiffraction_CreateGausSource( ObjectData ,&
             BeamWaistDiam,&
			 Nx,&
			 Ny,&
			 dx,&
			 dy,&
			 wavelength,&
			 BeamID,&
			 status) result (hresult)
        
        use Diffraction_Types

        implicit none
        
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        
        REAL(8), intent(in) :: BeamWaistDiam ! 1/e^2 intensity diameter
        INTEGER(4), intent(in) :: Nx
        INTEGER(4), intent(in) :: Ny
        REAL(8), intent(in) :: dx
        REAL(8), intent(in) :: dy
        REAL(8), intent(in) :: wavelength
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            call GenGaussBeam(ObjectData%Wbeams(BeamID), &
                Nx, Ny, dx, dy, wavelength, BeamWaistDiam, power=1.0_pr)
            status = INTMETCOMNOERROR
        ENDIF
        
        hresult = S_OK
    end function

    ! IDiffraction_CreateTophatSource
    function IDiffraction_CreateTophatSource( ObjectData ,&
             BeamDiam,&
			 Nx,&
			 Ny,&
			 dx,&
			 dy,&
			 wavelength,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none

        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData

        REAL(8), intent(in) :: BeamDiam
        INTEGER(4), intent(in) :: Nx
        INTEGER(4), intent(in) :: Ny
        REAL(8), intent(in) :: dx
        REAL(8), intent(in) :: dy
        REAL(8), intent(in) :: wavelength
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult

        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            call GenTopHatBeam(ObjectData%Wbeams(BeamID), &
                Nx, Ny, dx, dy, wavelength, BeamDiam, power=1.0_pr)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ClearWavefront
    function IDiffraction_ClearWavefront( ObjectData ,&
             BeamID,&
			 status) result (hresult)
        use Diffraction_Types

        implicit none

        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call ClearWavefront(ObjectData%Wbeams(BeamID))
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMask
    function IDiffraction_ApplyMask( ObjectData ,&
             length_r,&
			 length_t,&
			 offset,&
			 direction,&
			 BeamID,&
			 status) result (hresult)
        
        use Diffraction_Types
        use IntMetRoutines

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: length_r
        REAL(8), intent(in) :: length_t
        REAL(8), intent(in) :: offset
        CHARACTER(*), intent(in) :: direction
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
        
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call ApplyMask(ObjectData%Wbeams(BeamID), length_r, length_t, offset, direction)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMaskMisalign
    function IDiffraction_ApplyMaskMisalign( ObjectData ,&
             length_r,&
			 length_t,&
			 offset,&
			 xc,&
			 yc,&
			 direction,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines
        use Kinds
        use Wavefronts

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: length_r
        REAL(8), intent(in) :: length_t
        REAL(8), intent(in) :: offset
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        CHARACTER(*), intent(in) :: direction
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call ApplyMask(ObjectData%Wbeams(BeamID), length_r, length_t, offset, direction, &
                x_align=xc, y_align=yc )
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMaskGeneral
    !! Subroutine for applying a more general mask consisting of two
    !! rectangles symmetric about the origin
   function IDiffraction_ApplyMaskGeneral( ObjectData ,&
             len_x,&
			 len_y,&
			 off_x,&
			 off_y,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines
        use Kinds
        use Wavefronts

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: len_x
        REAL(8), intent(in) :: len_y
        REAL(8), intent(in) :: off_x
        REAL(8), intent(in) :: off_y
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! call ApplyMaskGeneral
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call ApplyMaskGeneral(ObjectData%Wbeams(BeamID), len_x, len_y, off_x, off_y)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMaskPoly
    function IDiffraction_ApplyMaskPoly( ObjectData ,&
             vertices,&
             xc,&
			 yc,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: vertices
        DIMENSION vertices(1:,1:)
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
        
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call ApplyMaskPoly(ObjectData%Wbeams(BeamID), vertices, xc, yc)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMaskX
    !! Mask consisting of four squares set on an 'X' pattern
    !! radius = distance from origin to center of each of the four squares
    !! length = side of each square
    function IDiffraction_ApplyMaskX( ObjectData ,&
             radius,&
			 length,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: radius
        REAL(8), intent(in) :: length
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
        
        ! call ApplyMaskX
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call ApplyMaskX(ObjectData%Wbeams(BeamID), radius, length)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMaskXRounded
    function IDiffraction_ApplyMaskXRounded( ObjectData ,&
             radius,&
			 length,&
			 corner_rad,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: radius
        REAL(8), intent(in) :: length
        REAL(8), intent(in) :: corner_rad
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! call ApplyMask
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call ApplyMaskXRounded(ObjectData%Wbeams(BeamID), radius, length, corner_rad)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMaskRotate
    function IDiffraction_ApplyMaskRotate( ObjectData ,&
             length_r,&
			 length_t,&
			 offset,&
			 angle,&
			 xc,&
			 yc,&
			 direction,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: length_r
        REAL(8), intent(in) :: length_t
        REAL(8), intent(in) :: offset
        REAL(8), intent(in) :: angle
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        CHARACTER(*), intent(in) :: direction
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
        call ApplyMaskRotate(ObjectData%Wbeams(BeamID), &
            length_r, length_t, offset, angle, direction, xc, yc)

            status = INTMETCOMNOERROR

        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ApplyMaskTilt
    function IDiffraction_ApplyMaskTilt( ObjectData ,&
             length_r,&
			 length_t,&
			 offset,&
			 rotangle,&
			 xc,&
			 yc,&
			 tiltangle,&
			 tiltorientation,&
			 direction,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: length_r
        REAL(8), intent(in) :: length_t
        REAL(8), intent(in) :: offset
        REAL(8), intent(in) :: rotangle
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        REAL(8), intent(in) :: tiltangle
        REAL(8), intent(in) :: tiltorientation
        CHARACTER(*), intent(in) :: direction
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
        call ApplyMaskRotate(ObjectData%Wbeams(BeamID), &
            length_r, length_t, offset, rotangle, direction, xc, yc, &
            tiltAngle=tiltangle, tiltOrientation=tiltorientation)

            status = INTMETCOMNOERROR

        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_Propagate
    function IDiffraction_Propagate( ObjectData ,&
             Distance,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: Distance
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
        ! local variables
        REAL(pr) :: dxout, dyout

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            ! fix dxout and dyout to be same as current sampling
            call GetWavefrontSampling(ObjectData%Wbeams(BeamID), dx=dxout, dy=dyout)

            call Propagate(ObjectData%Wbeams(BeamID), Distance, &
                dxout, dyout, applyCurv=.TRUE. )
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_PropagateExt
    function IDiffraction_PropagateExt( ObjectData ,&
             Distance,&
			 dxout,&
			 dyout,&
			 applyCurv,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none

        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: Distance
        REAL(8), intent(in) :: dxout
        REAL(8), intent(in) :: dyout
        LOGICAL(4), intent(in) :: applyCurv
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            call Propagate(ObjectData%Wbeams(BeamID), Distance, &
                dxout, dyout, applyCurv)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_PropagateDefDx
    ! uses default dx and dy values
    function IDiffraction_PropagateDefDx( ObjectData ,&
             Distance,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: Distance
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! dx and dy will be changed to the default values:
        ! far field: z*wavelength/Gx
        ! near field: (1 + z*curvature)*dx
        ! default value of applyCurv = .FALSE.
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            call Propagate(ObjectData%Wbeams(BeamID), Distance)
            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_CornerCube
    ! create a CornerCube object and reflect wavefront on it
    function IDiffraction_CornerCube( ObjectData ,&
             ccsize,&
			 shape,&
			 xc,&
			 yc,&
			 spin,&
			 gapwidth,&
			 rotmatrix,&
             dihedral,&
			 edgelength,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use Corner_Cubes
        use Array_Routines
        use Error_Exit

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: ccsize
        CHARACTER(*), intent(in) :: shape
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        REAL(8), intent(in) :: spin
        REAL(8), intent(in) :: gapwidth
        DIMENSION gapwidth(1:,1:)
        REAL(8), intent(in) :: rotmatrix
        DIMENSION rotmatrix(1:,1:)
        REAL(8), intent(in) :: dihedral
        DIMENSION dihedral(1:,1:)
        REAL(8), intent(in) :: edgelength
        DIMENSION edgelength(1:,1:)
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        !! lcoal variables
        type(CornerCube)         :: ccube
        REAL(pr), DIMENSION(3,3) :: Rmatrix
        REAL(pr)                 :: gw1, gw2, gw3, dh1, dh2, dh3, el1, el2, el3

        ! initialize to no error
        status = INTMETCOMNOERROR

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ENDIF

        IF ( status .EQ. INTMETCOMNOERROR ) THEN
            ! parse input gapwidth
            IF ( SIZE(gapwidth,2).EQ.3 ) THEN
                gw1 = gapwidth(1,1)
                gw2 = gapwidth(1,2)
                gw3 = gapwidth(1,3)
            ELSEIF ( SIZE(gapwidth,1).EQ.3 ) THEN
                gw1 = gapwidth(1,1)
                gw2 = gapwidth(2,1)
                gw3 = gapwidth(3,1)
            ELSE
                status = INVALIDGAPWIDTH
            ENDIF
        ENDIF

        IF ( status .EQ. INTMETCOMNOERROR ) THEN
            ! parse input dihedral
            IF ( SIZE(dihedral,2).EQ.3 ) THEN
                dh1 = dihedral(1,1)
                dh2 = dihedral(1,2)
                dh3 = dihedral(1,3)
            ELSEIF ( SIZE(dihedral,1).EQ.3 ) THEN
                dh1 = dihedral(1,1)
                dh2 = dihedral(2,1)
                dh3 = dihedral(3,1)
            ELSE
                status = INVALIDDIHEDRAL
            ENDIF
        ENDIF

        IF ( status .EQ. INTMETCOMNOERROR ) THEN
            ! parse input edgelength
            IF ( SIZE(edgelength,2).EQ.3 ) THEN
                el1 = edgelength(1,1)
                el2 = edgelength(1,2)
                el3 = edgelength(1,3)
            ELSEIF ( SIZE(edgelength,1).EQ.3 ) THEN
                el1 = edgelength(1,1)
                el2 = edgelength(2,1)
                el3 = edgelength(3,1)
            ELSE
                status = INVALIDEDGELENGTH
            ENDIF
        ENDIF

        IF ( status .EQ. INTMETCOMNOERROR ) THEN
            ! check that input rotmatrix is 3 x 3
            IF ( SIZE(rotmatrix,1).NE.3 .OR. SIZE(rotmatrix,2).NE.3 ) THEN
                status = INVALIDROTMATRIX
            ENDIF
        ENDIF

        IF ( status .EQ. INTMETCOMNOERROR ) THEN
            ! create corner cube object and apply to wavefront
            call NewCornerCube( ccube,                  &
                ccsize, shape, xc, yc, spin, rotmatrix, &
                dh1, dh2, dh3,                          & ! dihedral
                gw1, gw2, gw3,                          & ! gap widths
                el1, el2, el3 )                           ! edge lengths

            call ReflectFromCornerCube( ccube, ObjectData%Wbeams(BeamID) )
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_focuslens
    ! result = focuslens(real(beam),imag(beam))
    ! FocusLens:
    !     propagate one focal length
    !     apply lens
    !     propagate to focal plane
    ! see also PropagateToFocalPlane
    function IDiffraction_focuslens( ObjectData, &
             focuslens_f,&
			 focuslens_D,&
			 dxout,&
			 dyout,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use IntMetRoutines

        implicit none

        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: focuslens_f
        REAL(8), intent(in) :: focuslens_D
        REAL(8), intent(in) :: dxout
        REAL(8), intent(in) :: dyout
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
        
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            !! call FocusLens
            call FocusLens( ObjectData%Wbeams(BeamID), &
                focuslens_f, focuslens_D, dxout, dyout )

            status = INTMETCOMNOERROR

        ENDIF
                
        hresult = S_OK
    end function

    ! IDiffraction_ThinLens
    function IDiffraction_ThinLens( ObjectData ,&
             focallength,&
			 diameter,&
			 xc,&
			 yc,&
			 BeamID,&
			 status) result (hresult)
        
        use Diffraction_Types
        use Optics_Routines

        implicit none
        
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: focallength
        REAL(8), intent(in) :: diameter
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! hard wire for now
        LOGICAL applycurv
        applycurv = .FALSE.

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            !! Apply thin lens to wavefront
            call ThinLens( ObjectData%Wbeams(BeamID), &
                focallength, diameter, xc, yc, applycurv )
            status = INTMETCOMNOERROR

        ENDIF        
        
        hresult = S_OK
    end function

    ! IDiffraction_ReflectAsphere
    function IDiffraction_ReflectAsphere( ObjectData ,&
             rCurv,&
			 conicConst,&
			 diam,&
			 xDecenter,&
			 yDecenter,&
			 incidenceAngle,&
			 azimuth,&
			 applyCurv,&
			 BeamID,&
			 status) result (hresult)
        
        use Diffraction_Types
        use Optics_Routines

        implicit none
        
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: rCurv
        REAL(8), intent(in) :: conicConst
        REAL(8), intent(in) :: diam
        REAL(8), intent(in) :: xDecenter
        REAL(8), intent(in) :: yDecenter
        REAL(8), intent(in) :: incidenceAngle
        REAL(8), intent(in) :: azimuth
        INTEGER(4), intent(in) :: applyCurv
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
        
        ! local variables
        type(asphere)   mirror
        logical         logAC

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            !! make asphere object from input values
            call SetAsphereProperties( mirror, rCurv, conicConst, diam, &
                                        xDecenter, yDecenter )

            IF (applyCurv .NE. 0) THEN
                logAC = .TRUE.
            ELSE
                logAC = .FALSE.
            ENDIF

            !! apply the mirror to the wavefront
            call Reflect( mirror, ObjectData%Wbeams(BeamID), &
                          incidenceAngle, azimuth, logAC )
        
            status = INTMETCOMNOERROR

        ENDIF        
 
        hresult = S_OK
    end function


    ! IDiffraction_PistonTiltFocus
    function IDiffraction_PistonTiltFocus( ObjectData ,&
             piston,&
			 xtilt,&
			 ytilt,&
			 focus,&
			 xc,&
			 yc,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none

        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: piston
        REAL(8), intent(in) :: xtilt
        REAL(8), intent(in) :: ytilt
        REAL(8), intent(in) :: focus
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID
        ELSE
            call PistonTiltFocus( ObjectData%Wbeams(BeamID), &
                 piston, xtilt, ytilt, focus, xc, yc )
        
            status = INTMETCOMNOERROR
        ENDIF        

        hresult = S_OK
    end function

    ! IDiffraction_PropagateToFocalPlane
    ! applylens, then propaget to focal plane
    function IDiffraction_PropagateToFocalPlane( ObjectData ,&
             lens_f,&
			 lens_D,&
			 dxout,&
			 dyout,&
			 applycurv,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        USE Optics_Routines

        implicit none

        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: lens_f
        REAL(8), intent(in) :: lens_D
        REAL(8), intent(in) :: dxout
        REAL(8), intent(in) :: dyout
        LOGICAL(4), intent(in) :: applycurv
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        Type(simpleLens)  :: focuslens

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            !! call PropagateToFocalPlane
            call SetLensProperties( focuslens, lens_f, lens_D )
            call PropagateToFocalPlane( focuslens, &
                ObjectData%Wbeams(BeamID), dxout, dyout, applyCurv=applycurv )

            status = INTMETCOMNOERROR

        ENDIF


        hresult = S_OK
    end function

    ! IDiffraction_ClipCirc
    function IDiffraction_ClipCirc( ObjectData ,&
             Diameter,&
			 xc,&
			 yc,&
			 obscDiam,&
			 BeamID,&
			 status) result (hresult)
    ! clip the wavefront by a circular aperture with
    ! Diameter = Diameter
    ! xc, yc = center coordinate of aperture
    ! obscDiam = Diameter of center obscuration for creating an annulus
    !            If obscDiam <= 0, then no center obscuration
    ! if obscDiam > Diameter, then apply only obscuration to wavefront

        use Diffraction_Types

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: Diameter
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        REAL(8), intent(in) :: obscDiam
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            IF ( obscDiam >= Diameter ) THEN
                ! apply only obscuration
                call ClipCirc( ObjectData%Wbeams(BeamID), obscDiam, xc, yc, obsc=.TRUE. )
            ELSEIF ( obscDiam <= 0 ) THEN
                ! apply only circular clipping
                call ClipCirc( ObjectData%Wbeams(BeamID), Diameter, xc, yc )
            ELSE
                ! apply both circ clip and obscuration to form annulus
                call ClipCirc( ObjectData%Wbeams(BeamID), Diameter, xc, yc )
                call ClipCirc( ObjectData%Wbeams(BeamID), obscDiam, xc, yc, obsc=.TRUE. )
            ENDIF

            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ClipRect
    function IDiffraction_ClipRect( ObjectData ,&
             len_x,&
			 len_y,&
			 off_x,&
			 off_y,&
			 obscFlag,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: len_x
        REAL(8), intent(in) :: len_y
        REAL(8), intent(in) :: off_x
        REAL(8), intent(in) :: off_y
        INTEGER(4), intent(in) :: obscFlag
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status 
        integer(LONG) hresult
        ! TODO:  Add implementation
 
        LOGICAL obsctmp

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            ! non-zero obscFlag => obscuration mask
            IF ( obscFlag .NE. 0 ) THEN
                obsctmp = .TRUE.
            ELSE 
                obsctmp = .FALSE.
            ENDIF

            call ClipRect( ObjectData%Wbeams(BeamID), len_x, len_y, off_x, off_y, obsc=obsctmp )

            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ClipPolygon
    ! vertices = 2 rows by n columns, 
    !   first row = x coordinate
    !   second row = y coordinate
    !   the list of vertices must be counterclockwise
    function IDiffraction_ClipPolygon( ObjectData ,&
             vertices,&
             xc,&
             yc,&
			 angle,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types
        use Polygons

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: vertices
        DIMENSION vertices(1:,1:)
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        REAL(8), intent(in) :: angle
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation

        type(polygon) tmppoly

        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            ! create the polygon and the wavefront
            tmppoly = vertices
            call ClipPoly( ObjectData%Wbeams(BeamID), tmppoly, xc, yc, angle)

            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ClipBWindow
    function IDiffraction_ClipBWindow( ObjectData ,&
             winT,&
			 winA,&
			 xc,&
			 yc,&
			 rotAngle,&
			 BeamID,&
			 status) result (hresult)

        use Diffraction_Types

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: winT
        REAL(8), intent(in) :: winA
        REAL(8), intent(in) :: xc
        REAL(8), intent(in) :: yc
        REAL(8), intent(in) :: rotAngle
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult

        ! TODO:  Add implementation
        IF ( BeamID < 0 .or. BeamID > MAX_NR_WFRONTS ) THEN
            status = INVALIDWFRONTID

        ELSE
            call ClipWindow( ObjectData%Wbeams(BeamID), winT, winA, xc, yc, rotAngle )

            status = INTMETCOMNOERROR
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_AddThreadCommand
    function IDiffraction_AddThreadCommand( ObjectData ,&
             ThreadID,&
			 CmdID,&
			 rP,&
			 iP,&
			 pwAID,&
			 pwBID,&
			 strP,&
			 flagStatus) result (hresult)
        
        use Diffraction_Types
        use ThreadWrappers
        use Wavefronts
        
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(in) :: CmdID
        REAL(8), intent(in) :: rP
        DIMENSION rP(1:,1:)
        INTEGER(4), intent(in) :: iP
        DIMENSION iP(1:,1:)
        INTEGER(4), intent(in) :: pwAID
        INTEGER(4), intent(in) :: pwBID
        CHARACTER(*), intent(in) :: strP
        INTEGER(4), intent(out) :: flagStatus
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! ThreadID is index to ExecuteCommandsThreadParameters Tcmds
        ! CmdID is a command definition CMD_... defined in threadwrappers.f90
        ! rP is a vector of real valued parameters, specific for CmdID
        ! iP ...
        ! rP and iP are 2-d for Matlab compliance
        ! pwAID is index to wavefront array Wbeams()
        ! pwBID ...
        ! strP is string parameter, specific for CmdID
        ! flagStatus is the return value error code 0 = no error
        
        type(Wavefront), POINTER  :: pwA, pwB

        ! initialize to no error
        flagStatus = STATUS_NOERROR

        ! test valid inputs
        IF (ThreadID <= 0 .OR. ThreadID > MAX_NR_THREADS) THEN
            flagStatus = INVALIDTHREADID
        ENDIF

        ! assign temporary pointers to the wavefronts, or NULL
        IF (pwAID < 0 .OR. pwAID > MAX_NR_WFRONTS) THEN
            flagStatus = INVALIDWFRONTID
        ELSEIF (pwAID .EQ. 0) THEN
            NULLIFY(pwA)
        ELSE
            pwA => ObjectData%Wbeams(pwAID)
        ENDIF
        IF (pwBID < 0 .OR. pwBID > MAX_NR_WFRONTS) THEN
            flagStatus = INVALIDWFRONTID
        ELSEIF (pwBID .EQ. 0) THEN
            NULLIFY(pwB)
        ELSE
            pwB => ObjectData%Wbeams(pwBID)
        ENDIF

        ! if okay, add a new command to the linked list Tcmds(ThreadID)
        IF (flagStatus == STATUS_NOERROR) THEN
            flagStatus = NewCommand(ObjectData%Tcmds(ThreadID)%firstCommand, &
                            CmdID, rP, iP, pwA, pwB, strP)

        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ExecuteThread
    function IDiffraction_ExecuteThread( ObjectData ,&
             ThreadID,&
			 flagStatus) result (hresult)
        
        use Diffraction_Types
        use ThreadWrappers
        use DFWIN
        USE DFLIB

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(out) :: flagStatus
        integer(LONG) hresult
        ! TODO:  Add implementation

        INTEGER(4) :: hThread, thread_id
        flagStatus = STATUS_NOERROR

        ! test valid inputs
        IF (ThreadID <= 0 .OR. ThreadID > MAX_NR_THREADS) THEN
            flagStatus = INVALIDTHREADID
        ENDIF

        IF (flagStatus == STATUS_NOERROR) THEN

            ! call the thread
            hThread = CreateThread(                      &
                        0, 0,                            &
                        LOC(ExecuteCommandsThread),      &
                        LOC(ObjectData%Tcmds(ThreadID)), &
                        0,                               &
                        LOC(thread_id) ) 

            ObjectData%Tcmds(ThreadID)%hThread = hThread
            flagStatus = ObjectData%Tcmds(ThreadID)%flagStatus

        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_ThreadStatus
    function IDiffraction_ThreadStatus( ObjectData ,&
             ThreadID,&
			 flagStatus) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(out) :: flagStatus
        integer(LONG) hresult
        ! TODO:  Add implementation

        ! test valid inputs
        IF (ThreadID <= 0 .OR. ThreadID > MAX_NR_THREADS) THEN
            flagStatus = INVALIDTHREADID
        ELSE
            flagStatus = ObjectData%Tcmds(ThreadID)%flagStatus
        ENDIF

        hresult = S_OK
    end function

    ! IDiffraction_WaitForThread
    function IDiffraction_WaitForThread( ObjectData ,&
             ThreadID,&
			 TimeOutMS,&
			 WaitStatus) result (hresult)

        use Diffraction_Types
        use DFWIN
        USE DFLIB

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(in) :: TimeOutMS
        INTEGER(4), intent(out) :: WaitStatus
        integer(LONG) hresult
        ! TODO:  Add implementation

        WaitStatus = WaitForSingleObject(       &
            ObjectData%Tcmds(ThreadID)%hThread, TimeOutMS )

        hresult = S_OK
    end function

    ! IDiffraction_WaitForMultipleThreads
    function IDiffraction_WaitForMultipleThreads( ObjectData ,&
             pThreadID,&
			 TimeOutMS,&
			 flagWaitAll,&
			 WaitStatus) result (hresult)

        use Diffraction_Types
        use DFWIN
        USE DFLIB

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: pThreadID
        DIMENSION pThreadID(1:,1:)
        INTEGER(4), intent(in) :: TimeOutMS
        INTEGER(4), intent(in) :: flagWaitAll
        INTEGER(4), intent(out) :: WaitStatus
        integer(LONG) hresult
        ! TODO:  Add implementation

        INTEGER(HANDLE), DIMENSION(:), ALLOCATABLE :: hTList
        INTEGER(PHANDLE) :: lpHandles
        INTEGER                               :: nx, ny, i, j

        ny = size(pThreadID,1)
        nx = size(pThreadID,2)

        ALLOCATE(hTList(ny*nx))
        lpHandles = LOC(hTList)

        ! copy handles corresponding to ThreadID list into array
        do i = 1, ny
            do j = 1, nx
                hTList(i + (j-1)*ny) = ObjectData%Tcmds(pThreadID(i,j))%hThread
            enddo
        enddo

        WaitStatus = WaitForMultipleObjects(nx*ny, lpHandles, flagWaitAll, TimeOutMS)

        hresult = S_OK
    end function

    ! IDiffraction_ClearThread
    function IDiffraction_ClearThread( ObjectData ,&
             ThreadID,&
			 flagStatus) result (hresult)

        use Diffraction_Types
        use ThreadWrappers

        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(out) :: flagStatus
        integer(LONG) hresult
        ! TODO:  Add implementation

 
        ! test valid inputs
        IF (ThreadID < 0 .OR. ThreadID > MAX_NR_THREADS) THEN
            flagStatus = INVALIDTHREADID
        ELSEIF (ThreadID .EQ. 0) THEN
            ! single threaded methods are being used, do nothing
            flagStatus = STATUS_NOERROR
        ELSE
            call ClearCommandList(ObjectData%Tcmds(ThreadID)%firstCommand)
            NULLIFY(ObjectData%Tcmds(ThreadID)%firstCommand)
            ObjectData%Tcmds(ThreadID)%flagStatus = STATUS_NOERROR
            flagStatus = STATUS_NOERROR
        END IF

        hresult = S_OK
    end function

    ! IDiffraction_get_ThreadOPD
    function IDiffraction_get_ThreadOPD( ObjectData ,&
             ThreadID,&
			 OPDID,&
			 VALUE) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(in) :: OPDID
        REAL(8), intent(out) :: VALUE
        integer(LONG) hresult
        ! TODO:  Add implementation

        VALUE = 0.0_pr
        ! test valid inputs
        IF (ThreadID <= 0 .OR. ThreadID > MAX_NR_THREADS) THEN
            ObjectData%Tcmds(ThreadID)%flagStatus = INVALIDTHREADID
        ELSEIF (OPDID <= 0 .OR. OPDID > MAX_NR_OUTPUTS) THEN
            ObjectData%Tcmds(ThreadID)%flagStatus = INVALIDOUTPUTID
        ELSE
            VALUE = ObjectData%Tcmds(ThreadID)%outputOPD(OPDID)
        ENDIF

        hresult = S_OK
    end function
    ! IDiffraction_get_ThreadPOW
    function IDiffraction_get_ThreadPOW( ObjectData ,&
             ThreadID,&
			 POWID,&
			 VALUE) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(in) :: POWID
        REAL(8), intent(out) :: VALUE
        integer(LONG) hresult
        ! TODO:  Add implementation

        VALUE = 0.0_pr
        ! test valid inputs
        IF (ThreadID <= 0 .OR. ThreadID > MAX_NR_THREADS) THEN
            ObjectData%Tcmds(ThreadID)%flagStatus = INVALIDTHREADID
        ELSEIF (POWID <= 0 .OR. POWID > MAX_NR_OUTPUTS) THEN
            ObjectData%Tcmds(ThreadID)%flagStatus = INVALIDOUTPUTID
        ELSE
            VALUE = ObjectData%Tcmds(ThreadID)%outputPOW(POWID)
        ENDIF

        hresult = S_OK
    end function
