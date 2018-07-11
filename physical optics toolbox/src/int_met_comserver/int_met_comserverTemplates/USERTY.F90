!
!  @GLOBAL.SOURCE@ - This module contains user-defined class 
!  definitions and methods
!

module @INSTANCETYPE.USE@

    type @INSTANCETYPE.NAME@
        sequence
        !  TODO:  Add fields and remove "dummy"
        integer dummy
    end type

  contains

    !
    !  Constructor
    !
    function @INSTANCETYPE.CONSTRUCTOR@( ObjectData ) result (hresult)
        use dfwinty
        implicit none
        type(@INSTANCETYPE.NAME@) ObjectData
        !dec$ attributes reference :: ObjectData
        integer(LONG) hresult

        hresult = S_OK
        !  TODO:  Add field initialization code
    end function

    !
    !  Destructor
    !
    subroutine @INSTANCETYPE.DESTRUCTOR@( ObjectData )
        implicit none
        type(@INSTANCETYPE.NAME@) ObjectData
        !dec$ attributes reference :: ObjectData
        !  TODO:  Add field cleanup code
    end subroutine

end module

@($INSTANCENAME = CopyStr(INSTANCETYPE.NAME)@)
@#PER INTERFACE
@#template "userim.f90|U@INTERFACE.NAME@.f90|U@INTERFACE.NAME@+.f90"
@#ENDPER
