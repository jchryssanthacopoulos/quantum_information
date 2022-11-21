! Exercise 1
! ==========
!
! This program multiplies two matrices using the specified method and returns the time taken
!
! Command-line arguments:
!   mat_mul_method: Method to use to multiply the matrices (options are matmul, row-col, and col-row)
!     (default = `matmul`)
!   num_rows: Number of rows of the first matrix (default = 10)
!   num_cols: Number of columns of the second matrix (default = 10)
!   num_inner_dim: Number of inner dimensions (i.e., number of columns of first matrix, rows of second matrix)
!     (default = 10)
!   debug: Flag indicating whether to run in debug mode, which prints the inputs and outputs to the screen
!     (default = False)
!
! Raises:
!   Error if input dimensions are not positive integers
!
! Returns:
!   Time elapsed to multiply the matrices
!   If debug flag passed, input and output matrices are also displayed
!


! this module parses the command-line arguments needed to multiply the matrices
module arg_parse_mat_mul
    implicit none

    character(len=20) mat_mul_method
    character(len=20) char_dimensions(3)
    logical debug_mode

    ! default paramters
    character(len=50), parameter :: mat_mul_method_default = "matmul"
    character(len=50), parameter :: num_rows_default = "10"
    character(len=50), parameter :: num_cols_default = "10"
    character(len=50), parameter :: num_inner_dim_default = "10"
    logical, parameter :: debug_mode_default = .false.

contains
    ! parse command-line arguments
    subroutine parse_cmd_args
        implicit none

        integer ii
        integer num_args
        character(len=32) arg

        ! set defaults
        char_dimensions(1) = num_cols_default
        char_dimensions(2) = num_rows_default
        char_dimensions(3) = num_inner_dim_default
        mat_mul_method = mat_mul_method_default
        debug_mode = debug_mode_default

        num_args = command_argument_count()

        ii = 0

        do while(ii <= num_args)
            call get_command_argument(ii, arg)

            select case (arg)
                case ('--mat_mul_method')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) mat_mul_method
                        ii = ii + 1
                    end if

                case ('--num_rows')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) char_dimensions(1)
                        ii = ii + 1
                    end if

                case ('--num_cols')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) char_dimensions(2)
                        ii = ii + 1
                    end if

                case ('--num_inner_dim')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) char_dimensions(3)
                        ii = ii + 1
                    end if

                case ('-d', '--debug')
                    debug_mode = .true.
            end select

            ii = ii + 1
        end do

    end subroutine

end module


program exercise_1
    use arg_parse_mat_mul
    use arg_validate
    use mat_ops
    implicit none

    integer*4 dims(3)
    character(len=8) arg_char
    real*8 start, finish
    real*8, dimension(:, :), allocatable :: matrixA, matrixB, matprod

    ! parse command-line options
    call parse_cmd_args()

    ! check multiplication method is valid
    if (mat_mul_method .ne. "matmul" .and. mat_mul_method .ne. "row-col" .and. mat_mul_method .ne. "col-row") then
        print *, "Invalid matrix multiplication method. Options are: matmul, row-col, col-row"
        stop
    end if

    ! check if dimensions are integers
    call check_dims_integers(char_dimensions, dims, 3)
    if (status /= 0) then
        print *, "Dimensions need to be integers!"
        stop
    end if

    ! check if dimensions are positive
    call check_dims_positive(dims, 3)
    if (status /= 0) then
        print *, "Dimensions must be greater than zero!"
        stop
    end if

    ! print command-line options
    print *, "mat_mul_method = ", mat_mul_method
    write (arg_char, "(i8)") dims(1)
    print *, "num_rows = ", adjustl(arg_char)
    write (arg_char, "(i8)") dims(2)
    print *, "num_cols = ", adjustl(arg_char)
    write (arg_char, "(i8)") dims(3)
    print *, "num_inner_dim = ", adjustl(arg_char)

    ! allocate memory for matrices
    allocate(matrixA(dims(1), dims(3)))
    allocate(matrixB(dims(3), dims(2)))
    allocate(matprod(dims(1), dims(2)))

    ! initialize input matrices with random data
    call random_number(matrixA)
    call random_number(matrixB)

    if (debug_mode) then
        print *, "Running in debug mode ..."
        print *, "Matrix A = "
        call print_real_matrix(matrixA, dims(1), dims(3))
        print *, "Matrix B = "
        call print_real_matrix(matrixB, dims(3), dims(2))
    end if

    ! call intrinsic method
    call cpu_time(start)
        if (mat_mul_method .eq. "matmul") then
            matprod = matmul(matrixA, matrixB)
        else if (mat_mul_method .eq. "row-col") then
            call matmul_row_col(matrixA, matrixB, matprod, dims)
        else
            call matmul_col_row(matrixA, matrixB, matprod, dims)
        end if
    call cpu_time(finish)

    if (debug_mode) then
        print *, "Product ="
        call print_real_matrix(matprod, dims(1), dims(2))
    end if
    print "('Elapsed time = ', es16.10)", finish - start

    ! free memory
    deallocate(matrixA, matrixB, matprod)

end program
