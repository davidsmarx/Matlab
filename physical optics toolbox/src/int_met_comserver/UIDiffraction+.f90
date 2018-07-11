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
        hresult = S_OK
    end function
    ! IDiffraction_CopyWavefrontAmp
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
        hresult = S_OK
    end function
    ! IDiffraction_ReadWavefrontUNF
    function IDiffraction_ReadWavefrontUNF( ObjectData ,&
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
        REAL(8), intent(in) :: BeamWaistDiam
        INTEGER(4), intent(in) :: Nx
        INTEGER(4), intent(in) :: Ny
        REAL(8), intent(in) :: dx
        REAL(8), intent(in) :: dy
        REAL(8), intent(in) :: wavelength
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
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
        hresult = S_OK
    end function
    ! IDiffraction_ApplyMaskGeneral
    function IDiffraction_ApplyMaskGeneral( ObjectData ,&
             len_x,&
			 len_y,&
			 off_x,&
			 off_y,&
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
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
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
        hresult = S_OK
    end function
    ! IDiffraction_ApplyMaskX
    function IDiffraction_ApplyMaskX( ObjectData ,&
             radius,&
			 length,&
			 BeamID,&
			 status) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        REAL(8), intent(in) :: radius
        REAL(8), intent(in) :: length
        INTEGER(4), intent(in) :: BeamID
        INTEGER(4), intent(out) :: status
        integer(LONG) hresult
        ! TODO:  Add implementation
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
        hresult = S_OK
    end function
    ! IDiffraction_PropagateDefDx
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
        hresult = S_OK
    end function
    ! IDiffraction_CornerCube
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
        hresult = S_OK
    end function
    ! IDiffraction_focuslens
    function IDiffraction_focuslens( ObjectData ,&
             focuslens_f,&
			 focuslens_D,&
			 dxout,&
			 dyout,&
			 BeamID,&
			 status) result (hresult)
        use Diffraction_Types
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
        hresult = S_OK
    end function
    ! IDiffraction_PropagateToFocalPlane
    function IDiffraction_PropagateToFocalPlane( ObjectData ,&
             lens_f,&
			 lens_D,&
			 dxout,&
			 dyout,&
			 applycurv,&
			 BeamID,&
			 status) result (hresult)
        use Diffraction_Types
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
        hresult = S_OK
    end function
    ! IDiffraction_ClipPolygon
    function IDiffraction_ClipPolygon( ObjectData ,&
             vertices,&
			 xc,&
			 yc,&
			 angle,&
			 BeamID,&
			 status) result (hresult)
        use Diffraction_Types
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
        hresult = S_OK
    end function
    ! IDiffraction_ExecuteThread
    function IDiffraction_ExecuteThread( ObjectData ,&
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
        hresult = S_OK
    end function
    ! IDiffraction_WaitForThread
    function IDiffraction_WaitForThread( ObjectData ,&
             ThreadID,&
			 TimeOutMS,&
			 WaitStatus) result (hresult)
        use Diffraction_Types
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        INTEGER(4), intent(in) :: ThreadID
        INTEGER(4), intent(in) :: TimeOutMS
        INTEGER(4), intent(out) :: WaitStatus
        integer(LONG) hresult
        ! TODO:  Add implementation
        hresult = S_OK
    end function
    ! IDiffraction_WaitForMultipleThreads
    function IDiffraction_WaitForMultipleThreads( ObjectData ,&
             pThreadID,&
			 TimeOutMS,&
			 flagWaitAll,&
			 WaitStatus) result (hresult)
        use Diffraction_Types
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
        hresult = S_OK
    end function
    ! IDiffraction_ClearThread
    function IDiffraction_ClearThread( ObjectData ,&
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
        hresult = S_OK
    end function
