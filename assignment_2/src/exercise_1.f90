program exercise_1
    implicit none

    logical verbose
    integer i
    character(len=32) arg
    integer*4 a, b

    verbose = .false.
    a = 10
    b = 20

    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        select case (arg)
            case ('-v', '--verbose')
                verbose = .true.
        end select
    end do

    if (verbose) then
        call checkpoint(a, b)
    else
        print *, "Nothing to print. Exiting ..."
    end if

contains
    subroutine checkpoint(a, b)
        integer*4 a, b
        print *, "Entering checkpoint ..."
        print "('Value of variable a is ', i2)", a
        print "('Value of variable b is ', i2)", b
    end subroutine 

end program
