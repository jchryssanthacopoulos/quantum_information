! this module parses the command-line arguments needed to solve the Ising model in one dimension
module arg_parse
    implicit none

    integer ii

    integer N
    integer max_iter
    real*8 lambda
    real*8 thres
    character(len=50) diag_method
    logical debug

    ! default parameters
    integer, parameter :: N_default = 10
    integer, parameter :: max_iter_default = 50
    real*8, parameter :: lambda_default = 1.0
    real*8, parameter :: thres_default = 1d-10
    character(len=50), parameter :: diag_method_default = "zheev"
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
        max_iter = max_iter_default
        lambda = lambda_default
        thres = thres_default
        diag_method = diag_method_default
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

                case ("--max_iter")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) max_iter
                        ii = ii + 1
                    end if

                case ("--lambda")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) lambda
                        ii = ii + 1
                    end if

                case ("--thres")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) thres
                        ii = ii + 1
                    end if

                case ("--diag_method")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, "(A)") diag_method
                        ii = ii + 1
                    end if

                case ("-d", "--debug")
                    debug = .true.

            end select

            ii = ii + 1
        end do

    end subroutine

end module
