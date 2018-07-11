!
!  UDiffractionTY.f90 - This module contains user-defined class 
!  definitions and methods
!

module Diffraction_USE

    type Diffraction_InstanceData
        sequence
        !  TODO:  Add fields and remove "dummy"
        integer dummy
    end type

  contains

    !
    !  Constructor
    !
    function Diffraction_CONSTRUCTOR( ObjectData ) result (hresult)
        use dfwinty
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        integer(LONG) hresult

        hresult = S_OK
        !  TODO:  Add field initialization code
    end function

    !
    !  Destructor
    !
    subroutine Diffraction_DESTRUCTOR( ObjectData )
        implicit none
        type(Diffraction_InstanceData) ObjectData
        !dec$ attributes reference :: ObjectData
        !  TODO:  Add field cleanup code
    end subroutine

end module



