MODULE ThreadWrappers


    USE Kinds
    USE Wavefronts

    IMPLICIT NONE

    PUBLIC

    !! constant definitions for ExecuteCommandThreadParameters
    INTEGER, PARAMETER :: EXECUTESCRIPT_SUCCESS    = 0
    INTEGER, PARAMETER :: EXECUTESCRIPT_INPROGRESS = 1

    !! constant definitions for flagStatus
    INTEGER, PARAMETER :: STATUS_NOERROR    =  0
    INTEGER, PARAMETER :: STATUS_WORKING    =  1
    INTEGER, PARAMETER :: STATUS_EXITTHREAD = -1
    INTEGER, PARAMETER :: STATUS_ERROR      = -2
    INTEGER, PARAMETER :: STATUS_CMDIDERROR = -3

    !! command definitions
    INTEGER, PARAMETER :: CMD_WriteWavefrontUNF  = 1
    INTEGER, PARAMETER :: CMD_ReadWavefrontUNF   = 2
    INTEGER, PARAMETER :: CMD_ApplyMaskRotate    = 3
    INTEGER, PARAMETER :: CMD_CreateGausSource   = 4
    INTEGER, PARAMETER :: CMD_CreateTophatSource = 5
    INTEGER, PARAMETER :: CMD_Propagate          = 6
    INTEGER, PARAMETER :: CMD_CornerCube         = 7
    INTEGER, PARAMETER :: CMD_focuslens          = 8
    INTEGER, PARAMETER :: CMD_ClipCirc           = 9
    INTEGER, PARAMETER :: CMD_Heterodyne         = 10
    INTEGER, PARAMETER :: CMD_CopyWavefront      = 11
    INTEGER, PARAMETER :: CMD_WavefrontPower     = 12
    INTEGER, PARAMETER :: CMD_ClipPoly           = 13
    INTEGER, PARAMETER :: CMD_ApplyMaskPoly      = 14
    INTEGER, PARAMETER :: CMD_SumWavefront       = 15
    INTEGER, PARAMETER :: CMD_ThinLens           = 16
    INTEGER, PARAMETER :: CMD_ApplyMaskTilt      = 17
    INTEGER, PARAMETER :: CMD_CopyWavefrontAmp   = 18


    TYPE :: CommandStep
    ! linked list of commands, with members flexible for any sort
    ! of list of parameters. When opening a new thread, pass the
    ! address of the first command in the linked list and a StatusFlag
    ! When the StatusFlag
    ! is set, the thread will step through the list and execute each
    ! command step, then reset the StatusFlag when complete.
        sequence
        INTEGER*4                              :: iCommand
        REAL(pr), DIMENSION(:,:), POINTER      :: rP ! real parameters
        INTEGER*4, DIMENSION(:,:), POINTER     :: iP ! integer parameters
        TYPE(Wavefront), POINTER               :: pwA, pwB
        CHARACTER(256)                         :: strP ! string parameteres
        TYPE(CommandStep), POINTER             :: next
    END TYPE CommandStep

    TYPE :: ExecuteCommandsThreadParameters
        sequence
        TYPE(CommandStep), POINTER      :: firstCommand
        INTEGER(4)                      :: hThread   ! thread handle returned by CreateThread
        INTEGER                         :: flagStatus
        REAL(pr), DIMENSION(:), POINTER :: outputOPD ! vectors for storing calculation results
        REAL(pr), DIMENSION(:), POINTER :: outputPOW
    END TYPE ExecuteCommandsThreadParameters


CONTAINS
    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!
    ! subroutine to allocate a command and add it to the end of the list
    FUNCTION NewCommand (pCmd, iCmd, rP, iP, pwA, pwB, strP) result (status)

    IMPLICIT NONE

    type(CommandStep), POINTER                :: pCmd ! firstCommand in the linked list
    integer, INTENT(IN)                       :: iCmd
    real(pr), DIMENSION(:,:), INTENT(IN)      :: rP ! 2-d for Matlab compliance
    integer*4, DIMENSION(:,:), INTENT(IN)     :: iP ! 2-d for matlab compliance
    type(Wavefront), POINTER                  :: pwA, pwB
    character(*), INTENT(IN)                  :: strP
    integer                                   :: status
    ! status returns the number of commands in the linked list, or <0 if error

    type(CommandStep), pointer :: new_cmd, lastcmd
    integer ccnt

    status = STATUS_NOERROR ! intitialize, change the status if a problem occurs

    ! allocate a new command and load values
    allocate(new_cmd)
    new_cmd%iCommand = iCmd

    ! copy real and integer data
    allocate(new_cmd%rP(size(rP,1),size(rP,2)))
    new_cmd%rP = rP
    allocate(new_cmd%iP(size(iP,1),size(iP,2)))
    new_cmd%iP = iP 

    new_cmd%pwA => pwA
    new_cmd%pwB => pwB

    new_cmd%strP = " "
    IF ( len(strP) > len(new_cmd%strP) ) THEN
        status = STATUS_ERROR
    ELSE
        new_cmd%strP(1:len(strP)) = strP(1:len(strP))
    END IF

    NULLIFY(new_cmd%next)

    ! cycle down the linked list to the last command, then add the
    ! new Cmd by setting the last command%next to the new command    
    IF ( status == STATUS_NOERROR ) THEN

        IF ( .not. associated(pCmd) ) THEN
            ! this must be the first command in the list
            pCmd => new_cmd
            ccnt = 1
        ELSE
            ! go down the list to the last command, and add the new_cmd
            lastcmd => pCmd
            ccnt = 1
            DO WHILE ( associated(lastcmd%next) )
                lastcmd => lastcmd%next
                ccnt = ccnt + 1
            END DO
            ! assign the new command to the parent
            lastcmd%next => new_cmd
            ccnt = ccnt + 1
        END IF ! associated

    ELSE ! there was an error, clear new_cmd
        call ClearCommand(new_cmd)

    END IF ! no error

    END FUNCTION NewCommand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! clear a single command
    SUBROUTINE ClearCommand (pCmd)

    IMPLICIT NONE

    type(CommandStep), POINTER  :: pCmd

    deallocate (pCmd%rP, pCmd%iP)
    deallocate (pCmd)

    END SUBROUTINE ClearCommand

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! clear a whole list of commands
    SUBROUTINE ClearCommandList (pCmd)

    IMPLICIT NONE

    type(CommandStep), POINTER  :: pCmd
    type(CommandStep), POINTER  :: thiscmd, nextcmd

    ! store next command and deallocate current
    nextcmd => pCmd
    DO WHILE ( associated(nextcmd) )
        thiscmd => nextcmd
        nextcmd => thiscmd%next
        deallocate(thiscmd%rP, thiscmd%iP)
        deallocate(thiscmd)
    END DO

    END SUBROUTINE ClearCommandList


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
    !! Subroutine to execute the commands described in a linked list
    !! to use the thread, create a linked list of commands
    !! associate the first command in the list to the firstCommand
    !! set thread flagStatus = STATUS_WORKING
    !! call CreateThread
    SUBROUTINE ExecuteCommandsThread (ptp)
    
    USE dfwin
    
    IMPLICIT NONE

    TYPE(ExecuteCommandsThreadParameters), POINTER :: ptp
!DEC$ ATTRIBUTES VALUE :: ptp

    TYPE(CommandStep), POINTER :: cmdCurrent

    !!!
    cmdCurrent => ptp%firstCommand

    !!!! main thread loop, keeps running until error or told to stop
    !! flag values and what they mean:
    !        STATUS_NOERROR     nothing pending, wait
    !        STATUS_WORKING     execute next command
    !        STATUS_EXITTHREAD  finish executing commands, then exit thread
    !        STATUS_ERROR       error executing command, exit thread now
    

    IF (associated(cmdCurrent)) THEN

        ptp%flagStatus = STATUS_WORKING
        ! command loop
        DO WHILE (ptp%flagStatus.eq.STATUS_WORKING)

            ! execute current command
            call ExecuteCommand ! cmdCurrent

            ! flush the wavefronts for use by other threads and next command
!$OMP FLUSH (cmdCurrent%pwA, cmdCurrent%pwB)

            ! if no error, increment to next command and check if end of list
            if (ptp%flagStatus .eq. STATUS_WORKING) then
                
                if ( associated(cmdCurrent%next) ) then
                    ! increment cmdCurrent
                    cmdCurrent => cmdCurrent%next
                else
                    ! reached end of list with no error
                    ! set flag to signal thread is finished
                    ptp%flagStatus = STATUS_NOERROR
                end if

            end if ! if no error, otherwise loop with exit

        END DO ! command loop

    END IF ! valid first command

    ! exit the thread
    call ExitThread(ptp%flagStatus)


    CONTAINS
    
        !! subroutine the parses the command and calls the appropriate
        !! diffraction code
        !!INTEGER*4                           :: iCommand
        !!REAL(pr), DIMENSION(:), ALLOCATABLE :: rP
        !!TYPE(Wavefront), POINTER            :: pwA, pwB
        !!CHARACTER*64                        :: strP
        !!TYPE(CommandStep), POINTER          :: next
        SUBROUTINE ExecuteCommand
            
           USE IntMetRoutines
           USE Corner_Cubes
           USE Polygons
           USE Optics_Routines

           IMPLICIT NONE
            
           ! local variables
           REAL(pr)                          :: dxout, dyout
           type(CornerCube)                  :: ccube
           REAL(pr), DIMENSION(:,:), POINTER :: rP ! tmp pointer for convenience
           INTEGER                           :: ii, jj ! tmp variables
           type(polygon)                     :: tmppoly
           LOGICAL                           :: applycurv


           SELECT CASE (cmdCurrent%iCommand)

           case(CMD_WriteWavefrontUNF)
              call WriteWavefrontUNF(cmdCurrent%pwA, cmdCurrent%strP )

           case(CMD_ReadWavefrontUNF)
              call ReadWavefrontUNF(cmdCurrent%pwA, cmdCurrent%strP )

           case(CMD_ApplyMaskRotate)
              !!call ApplyMaskRotate(ObjectData%wvfront, length_r, length_t, offset, &
              !!                    angle, direction, xc, yc)
              call ApplyMaskRotate(cmdCurrent%pwA, &
                  cmdCurrent%rP(1,1), cmdCurrent%rP(1,2), cmdCurrent%rP(1,3), &
                  cmdCurrent%rP(1,4), cmdCurrent%strP, cmdCurrent%rP(1,5), cmdCurrent%rP(1,6))

           case(CMD_ApplyMaskTilt)
               !!call ApplyMaskRotate with additional tilt parameters to tilt the mask
               !!with respect to the propagation axis
               ! rP(1,1:8) = [length_r, length_t, offset, rotangle (about the propagation axis
               !              xc, yc, tiltangle (AOI), tiltOrientation (angle between tilt axis
               !              and x-axis)
               ! strP = 'ns' or 'ew'
               call ApplyMaskRotate(cmdCurrent%pwA, &
                  cmdCurrent%rP(1,1), cmdCurrent%rP(1,2), cmdCurrent%rP(1,3), &
                  cmdCurrent%rP(1,4), cmdCurrent%strP, cmdCurrent%rP(1,5), cmdCurrent%rP(1,6), &
                  cmdCurrent%rP(1,7), cmdCurrent%rP(1,8))


           case(CMD_ApplyMaskPoly)
               !!call ApplyMaskPoly(ObjectData%wvfront, vertices, xc, yc)
               ! rP(1,1:n) = x-coordinates of polygon
               ! rP(2,1:n) = y-coordinates of polygon
               ! rP(3,1)   = xc
               ! rP(4,1)   = yc
               call ApplyMaskPoly(cmdCurrent%pwA, cmdCurrent%rP(1:2,:), cmdCurrent%rP(3,1), cmdCurrent%rP(4,1))


           case(CMD_CreateGausSource)
               !!call GenGaussBeam(ObjectData%wvfront, &
                   !!Nx, Ny, dx, dy, wavelength, BeamWaistDiam)
               call GenGaussBeam(cmdCurrent%pwA, cmdCurrent%iP(1,1), cmdCurrent%iP(1,2), &
                   cmdCurrent%rP(1,1),cmdCurrent%rP(1,2),cmdCurrent%rP(1,3),cmdCurrent%rP(1,4),&
                   power=1.0_pr)
 
              
           case(CMD_CreateTophatSource)
               !!call GenTopHatBeam(ObjectData%wvfront, &
                    !!Nx, Ny, dx, dy, wavelength, BeamDiam)
               call GenTopHatBeam(cmdCurrent%pwA, cmdCurrent%iP(1,1), cmdCurrent%iP(1,2), &
                   cmdCurrent%rP(1,1),cmdCurrent%rP(1,2),cmdCurrent%rP(1,3),cmdCurrent%rP(1,4), &
                   power=1.0_pr)


           case(CMD_Propagate)
                ! rP(1,1) = distance
                ! rP(1,2:3) = [dxout dyout] (optional, default = current values)
                ! iP(1,1) = 1 => applyCurv = .TRUE., 0 => applyCurv = .FALSE.
                !
                rP => cmdCurrent%rP ! for ease of typing

                IF ( SIZE(rP,2) < 3 ) THEN
                     call GetWavefrontSampling(cmdCurrent%pwA, dx=dxout, dy=dyout)
                ELSE
                     dxout = rP(1,2)
                     dyout = rP(1,3)
                ENDIF

                applycurv = .FALSE. ! default
                IF ( cmdCurrent%iP(1,1) .eq. 1 ) applycurv = .TRUE.

                call Propagate(cmdCurrent%pwA, cmdCurrent%rP(1,1), &
                    dxout, dyout, applyCurv=applycurv )


           case(CMD_CornerCube)
                ! create corner cube object and apply to wavefront
                ! rP = [ccsize xc yc
                !       spin 0 0
                !       gapwidth(1:3)
                !       dihedral(1:3)
                !       edgelength(1:3)
                !       rotmatrix(1:3,1:3)]
                ! strP = shape
                ! apply to pwA
                                
                rP => cmdCurrent%rP ! just for ease of typing

                call NewCornerCube( ccube, rP(1,1), cmdCurrent%strP, rP(1,2), rP(1,3), rP(2,1), &
                                    rP(6:8,1:3), &
                                    rP(4,1), rP(4,2), rP(4,3), &
                                    rP(3,1), rP(3,2), rP(3,3), &
                                    rP(5,1), rP(5,2), rP(5,3) )

                call ReflectFromCornerCube( ccube, cmdCurrent%pwA )
           
           
           case(CMD_focuslens)

                !!call FocusLens( ObjectData%wvfront, focuslens_f, focuslens_D, dxout, dyout )
                call FocusLens( cmdCurrent%pwA, cmdCurrent%rP(1,1), cmdCurrent%rP(1,2),&
                                cmdCurrent%rP(1,3), cmdCurrent%rP(1,4) )


           case(CMD_ThinLens)
                ! apply thin lens to wavefront A
                ! rP(1,1) = focallength
                ! rP(1,2) = diameter
                ! rP(1,3) = xc
                ! rP(1,4) = yc
                ! applycurv is fixed = false
                call ThinLens( cmdCurrent%pwA, cmdCurrent%rP(1,1), cmdCurrent%rP(1,2), &
                               cmdCurrent%rP(1,3), cmdCurrent%rP(1,4), applyCurv=.FALSE. )


           case(CMD_ClipCirc)
               ! rP(1) = diameter, rP(2,3) = xc,yc, rP(4) = obscDiam
               IF ( cmdCurrent%rP(1,4) .GT. 0.0_pr ) THEN
                    IF ( cmdCurrent%rP(1,4) .GE. cmdCurrent%rP(1,1) ) THEN
                        call ClipCirc( cmdCurrent%pwA, cmdCurrent%rP(1,4), cmdCurrent%rP(1,2), cmdCurrent%rP(1,3), obsc=.TRUE. )
                    ELSE
                        call ClipCirc( cmdCurrent%pwA, cmdCurrent%rP(1,1), cmdCurrent%rP(1,2), cmdCurrent%rP(1,3) )
                        call ClipCirc( cmdCurrent%pwA, cmdCurrent%rP(1,4), cmdCurrent%rP(1,2), cmdCurrent%rP(1,3), obsc=.TRUE. )
                    ENDIF
               ELSE
                    call ClipCirc( cmdCurrent%pwA, cmdCurrent%rP(1,1), cmdCurrent%rP(1,2), cmdCurrent%rP(1,3) )
               ENDIF


           case(CMD_ClipPoly)
               ! rP(1,1:n) = x-coordinates of polygon
               ! rP(2,1:n) = y-coordinates of polygon
               ! rP(3,1)   = xc
               ! rP(4,1)   = yc
               ! rP(5,1)   = angle
               
               ! create the polygon and the wavefront
               tmppoly = cmdCurrent%rP(1:2,:)
               call ClipPoly( cmdCurrent%pwA, tmppoly, &
                    cmdCurrent%rP(3,1), cmdCurrent%rP(4,1), cmdCurrent%rP(5,1) )


           case(CMD_Heterodyne)
               ! heterodyne pwA and pwB
               ! rP(1) = xSize, rP(2) = ySize, rP(3,4) = (xc,yc)
               ! iP(1) = index to store opd result (outputOPD(iP(1)) = OPD)
               ! iP(2) = index to store opd result (outputPOW(iP(2)) = power)
               ! strP = detShape
               call Heterodyne( cmdCurrent%pwA, cmdCurrent%pwB, &
                    ptp%outputOPD(cmdCurrent%iP(1,1)), ptp%outputPOW(cmdCurrent%iP(1,2)), &
                    cmdCurrent%strP, &
                    cmdCurrent%rP(1,1), cmdCurrent%rP(1,2), cmdCurrent%rP(1,3), cmdCurrent%rP(1,4) )

           case (CMD_CopyWavefront)
               ! copy wavefront from pwA to pwB
               cmdCurrent%pwB = cmdCurrent%pwA

           case (CMD_CopyWavefrontAmp)
               ! copy the amplitude from pwA to pwB, so pwB will have its original
               ! phase with the amplitude of pwA
               call CopyWavefrontAmp(cmdCurrent%pwA, cmdCurrent%pwB)

           case (CMD_SumWavefront)
               ! cumulutive sum, a = a+b
               call WavefrontCumSum(cmdCurrent%pwA,cmdCurrent%pwB)

           case (CMD_WavefrontPower)
               ! output current wavefront power of
               ! pwA to:
               ! outputPOW(iP(1))
               ptp%outputPOW(cmdCurrent%iP(1,1)) = WavefrontPower(cmdCurrent%pwA)

           case default
              ptp%flagStatus = STATUS_CMDIDERROR

           END SELECT        


        END SUBROUTINE ExecuteCommand


    END SUBROUTINE ExecuteCommandsThread

END MODULE ThreadWrappers