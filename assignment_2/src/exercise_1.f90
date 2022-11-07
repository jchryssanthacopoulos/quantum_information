! Exercise 1
! ==========
!
! This program contains a subroutine to be used as a checkpoint for debugging, which can be entered by passing the
!   -d/--debug flag on the command-line
!
! Flags:
!   -d/--debug to enter debug checkpoint and inspect values of variables
!
! Returns:
!   Values of certain variables if debug flag passed in
!


program exercise_1
    implicit none

    logical debug_mode
    character(len=32) arg
    integer*4 ii, a, b

    debug_mode = .false.
    a = 10
    b = 20

    ! parse command-line options
    do ii = 1, command_argument_count()
        call get_command_argument(ii, arg)
        select case (arg)
            case ('-d', '--debug')
                debug_mode = .true.
        end select
    end do

    ! enter checkpoint if in debug mode
    if (debug_mode) then
        call checkpoint(a, b)
    else
        print *, "Nothing to print. Exiting ..."
    end if

contains
    ! debug checkpoint that prints values of certain variables
    !
    ! Inputs:
    !   a: First variable to print
    !   b: Second variable to print
    !
    subroutine checkpoint(a, b)
        integer*4 a, b
        print *, "Entering checkpoint ..."
        print "('Value of variable a is ', i2)", a
        print "('Value of variable b is ', i2)", b
    end subroutine 

end program
