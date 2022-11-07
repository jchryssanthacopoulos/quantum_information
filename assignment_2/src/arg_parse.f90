! this module contains utilities to parse command-line arguments and validate matrix inputs
module arg_parse
    implicit none

    integer*4 ii

    ! used for indicating status for argument parsing
    integer status

    ! command-line arguments
    logical debug_mode

contains
    ! parse command-line arguments
    subroutine parse_cmd_args
        implicit none

        character(len=32) arg

        do ii = 1, command_argument_count()
            call get_command_argument(ii, arg)
            select case (arg)
                case ('-d', '--debug')
                    debug_mode = .true.
            end select
        end do
    end subroutine

    ! checks whether the dimensions are integers
    !
    ! Inputs:
    !   char_input: Inputs to check
    !   dims: Where to save dimensions if integer conversion could be performed
    !   num_dims: Number of dimensions to check
    !
    subroutine check_dims_integers(char_input, dims, num_dims)
        implicit none

        integer*4 num_dims
        character(len=20) char_input(num_dims)
        integer*4 dims(num_dims)

        do ii = 1, num_dims
            call str2int(char_input(ii), dims(ii), status)
            if (status /= 0) then
                return
            end if
        end do
    end subroutine

    ! checks whether the dimensions are positive
    !
    ! Inputs:
    !   dims: Where to save dimensions if integer conversion could be performed
    !   num_dims: Number of dimensions to check
    !
    subroutine check_dims_positive(dims, num_dims)
        implicit none

        integer*4 num_dims
        integer*4 dims(num_dims)

        do ii = 1, num_dims
            if (dims(ii) < 1) then
                status = -1
                return
            end if
        end do
    end subroutine

    ! converts a string into an integer, returning a non-zero status code if conversion failed
    ! adapted from code found here:
    ! https://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90
    !
    ! Inputs:
    !   str: String to try to convert to integer
    !   int: Integer to save the result to
    !   status: Resulting status (if it is equal to zero, conversion was successful; non-zero otherwise)
    !
    subroutine str2int(str, int, status)
        implicit none

        character(len=*), intent(in) :: str
        integer, intent(out) :: int
        integer, intent(out) :: status

        ! read string into an integer, saving the resulting status
        read(str, *, iostat=status) int
    end subroutine

end module
