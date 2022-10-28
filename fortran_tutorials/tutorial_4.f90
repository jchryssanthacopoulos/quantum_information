program main
    implicit none

    logical verbose
    integer i
    character(len=32) arg

    verbose = .false.

    do i = 1, command_argument_count()
        call get_command_argument(i, arg)
        select case (arg)
            case ('-v', '--verbose')
                verbose = .true.
        end select
    end do

    if (verbose) then
        print *, "verbose on"
    else
        print *, "verbose off"
    end if

end program
