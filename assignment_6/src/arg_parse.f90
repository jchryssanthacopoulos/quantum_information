! this module parses the command-line arguments needed to generate a density matrix of a many-body quantum system
module arg_parse
    implicit none

    integer ii

    integer N, D, M
    character(len=50) system_type
    character(len=50) output_filename
    integer debug_level

    ! default parameters
    integer, parameter :: N_default = 2
    integer, parameter :: D_default = 2
    integer, parameter :: M_default = 1
    character(len=50), parameter :: system_type_default = "separable"
    character(len=50), parameter :: output_filename_default = "density.txt"
    integer, parameter :: debug_level_default = 0

    ! for displaying command-line arguments
    character(len=8) arg_char

    integer, parameter :: DEBUG_LEVEL_1 = 1
    integer, parameter :: DEBUG_LEVEL_2 = 2

contains
    ! parse command-line arguments
    subroutine parse_cmd_args
        implicit none

        integer num_args
        character(len=32) arg

        ! set defaults
        N = N_default
        D = D_default
        M = M_default
        system_type = system_type_default
        output_filename = output_filename_default
        debug_level = debug_level_default

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

                case ("--D")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) D
                        ii = ii + 1
                    end if

                case ("--M")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) M
                        ii = ii + 1
                    end if

                case ("--type")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, "(A)") system_type
                        ii = ii + 1
                    end if

                case ("--output_filename")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, "(A)") output_filename
                        ii = ii + 1
                    end if

                case ("--debug")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) debug_level
                        ii = ii + 1
                    end if

            end select

            ii = ii + 1
        end do

    end subroutine

end module
