@#IFANY METHOD
!
!  @GLOBAL.SOURCE@ - This file contains the implementation of the
!                    @INTERFACE.NAME@ methods
!
    @#PER METHOD
    ! @INTERFACE.NAME@_@METHOD.NAME@
    function @INTERFACE.NAME@_@METHOD.NAME@( ObjectData@#IFANY ARGUMENT ,@#ENDIF &
            @#PER ARGUMENT @ARGUMENT.NAME@@#PERPOSTFIX ",&
			"@#ENDPER ) result (hresult)
        use @CLASS.NAME@_Types
        implicit none
        type(@CLASS.$INSTANCENAME@) ObjectData
        !dec$ attributes reference :: ObjectData
        @#PER ARGUMENT
        @(= FortDefinition(ARGUMENT.DATATYPE, ARGUMENT.INTENT, ARGUMENT.ARRAY, ARGUMENT.OPTIONAL)@) :: @ARGUMENT.NAME@@[
        @(= FortDimension2(ARGUMENT.NAME, ARGUMENT.ARRAY, ARGUMENT.RANK, ARGUMENT.ASSUMEDSHAPE, ARGUMENT.LB1, ARGUMENT.UB1, ARGUMENT.LB2, ARGUMENT.UB2, ARGUMENT.LB3, ARGUMENT.UB3, ARGUMENT.LB4, ARGUMENT.UB4, ARGUMENT.LB5, ARGUMENT.UB5, ARGUMENT.LB6, ARGUMENT.UB6, ARGUMENT.LB7, ARGUMENT.UB7)@)@]
        @#ENDPER
        integer(LONG) hresult
        ! TODO:  Add implementation
        hresult = S_OK
    end function
    @#ENDPER
@#ENDIF
