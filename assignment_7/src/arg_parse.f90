! this module parses the command-line arguments needed to solve the Ising model in one dimension
module arg_parse
    implicit none

    integer ii

    integer N
    real*8 lambda
    character(len=50) output_filename
    logical debug

    ! default parameters
    integer, parameter :: N_default = 10
    integer, parameter :: lambda_default = 1.0
    character(len=50), parameter :: output_filename_default = "ising.txt"
    logical, parameter :: debug_default = .false.

    ! for displaying command-line arguments
    character(len=8) arg_char

contains
    ! parse command-line arguments
    subroutine parse_cmd_args
        implicit none

        integer num_args
        character(len=32) arg

        ! set defaults
        N = N_default
        lambda = lambda_default
        output_filename = output_filename_default
        debug = debug_default

        num_args = command_argument_count()

        ii = 0

        do while (ii <= num_args)
            call get_command_argument(ii, arg)

            select case (arg)
                case ("--N")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) N
                        ii = ii + 1
                    end if

                case ("--lambda")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) lambda
                        ii = ii + 1
                    end if

                case ("--output_filename")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, "(A)") output_filename
                        ii = ii + 1
                    end if

                case ("-d", "--debug")
                    debug = .true.

            end select

            ii = ii + 1
        end do

    end subroutine

end module
