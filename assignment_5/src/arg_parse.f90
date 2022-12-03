! this module parses the command-line arguments needed to solve the time-dependent Schrodinger equation
module arg_parse
    implicit none

    integer ii

    real*8 xmin, xmax, tmax
    integer num_x_pts, num_t_pts
    character(len=50) output_filename

    ! default parameters
    real*8, parameter :: xmin_default = -5.0
    real*8, parameter :: xmax_default = 5.0
    real*8, parameter :: tmax_default = 10.0
    integer, parameter :: num_x_pts_default = 1000
    integer, parameter :: num_t_pts_default = 100
    character(len=50), parameter :: output_filename_default = "solution.txt"

    ! for displaying command-line arguments
    character(len=8) arg_char

contains
    ! parse command-line arguments
    subroutine parse_cmd_args
        implicit none

        integer num_args
        character(len=32) arg

        ! set defaults
        xmin = xmin_default
        xmax = xmax_default
        tmax = tmax_default
        num_x_pts = num_x_pts_default
        num_t_pts = num_t_pts_default
        output_filename = output_filename_default

        num_args = command_argument_count()

        ii = 0

        do while (ii <= num_args)
            call get_command_argument(ii, arg)

            select case (arg)
                case ("--xmin")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) xmin
                        ii = ii + 1
                    end if

                case ("--xmax")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) xmax
                        ii = ii + 1
                    end if

                case ("--tmax")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) tmax
                        ii = ii + 1
                    end if

                case ("--num_x_pts")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) num_x_pts
                        ii = ii + 1
                    end if

                case ("--num_t_pts")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) num_t_pts
                        ii = ii + 1
                    end if

                case ("--output_filename")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, "(A)") output_filename
                        ii = ii + 1
                    end if

            end select

            ii = ii + 1
        end do

    end subroutine

end module
