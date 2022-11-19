! Exercise 1
! ==========
!
! This program implements two different ways of performing matrix multiplication and compares their results and
!   performance against the builtin method, matmul
!
! Inputs:
!   number of rows of first matrix
!   number of columns of second matrix
!   number of inner dimensions (i.e., number of columns of first matrix, rows of second matrix)
!
! Flags:
!   -d/--debug to display input and output matrices
!
! Raises:
!   Error if inputs are not positive integers
!
! Returns:
!   Elapsed time for each multiplication method and maximum absolute error
!   If debug flag set, input and output matrices are also displayed
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


! this module contains debug checkpoints to print input and output matrices
module debug
    implicit none

contains
    ! print input matrices to multiply
    !
    ! Inputs:
    !   A: First matrix to multiply
    !   B: Second matrix to multiply
    !   dims: Number of rows of A, columns of B, and inner dimension
    !
    subroutine print_input_matrices(A, B, dims)
        implicit none

        integer*4 dims(3)
        real*8 A(dims(1), dims(3)), B(dims(3), dims(2))

        print *, "Matrix A = "
        call print_matrix(A, dims(1), dims(3))
        print *, "Matrix B = "
        call print_matrix(B, dims(3), dims(2))
    end subroutine

    ! print matrix computed using given method
    !
    ! Inputs:
    !   A: Computed matrix
    !   nrows: Number of rows
    !   ncols: Number of columns
    !   method_name: Method used
    !
    subroutine print_matrix_for_method(A, nrows, ncols, method_name)
        implicit none

        integer*4 nrows, ncols
        real*8 A(nrows, ncols)
        character(len=*) method_name

        print *, "Product using ", method_name, " = "
        call print_matrix(A, nrows, ncols)
    end subroutine

    ! print matrix in a nice format
    !
    ! Inputs:
    !   M: Matrix to print
    !   nrows: Number of rows
    !   ncols: Number of columns
    !
    subroutine print_matrix(M, nrows, ncols)
        implicit none

        integer*4 ii, nrows, ncols
        real*8 M(nrows, ncols)

        do ii = 1, nrows
            print '(20f7.2)', M(ii, 1:ncols)
        end do
    end subroutine

end module


program exercise_1
    use arg_validate
    use arg_parse_mat_mul
    use debug
    use mat_ops
    implicit none

    integer*4 dims(3)
    real*8 start, finish
    real*8, dimension(:, :), allocatable :: matrixA, matrixB, matprod

    ! parse command-line options
    call parse_cmd_args()

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
    print *, "num_rows = ", dims(1)
    print *, "num_cols =", dims(2)
    print *, "num_inner_dim =", dims(3)

    ! allocate memory for matrices
    allocate(matrixA(dims(1), dims(3)))
    allocate(matrixB(dims(3), dims(2)))
    allocate(matprod(dims(1), dims(2)))

    ! initialize input matrices with random data
    call random_number(matrixA)
    call random_number(matrixB)

    if (debug_mode) then
        print *, "Running in debug mode ..."
        call print_input_matrices(matrixA, matrixB, dims)
    end if

    ! call intrinsic method
    call cpu_time(start)
        if (mat_mul_method .eq. "matmul") then
            matprod = matmul(matrixA, matrixB)
        else if (mat_mul_method .eq. "row-col") then
            call matmul_row_col(matrixA, matrixB, matprod, dims)
        else if (mat_mul_method .eq. "col-row") then
            call matmul_col_row(matrixA, matrixB, matprod, dims)
        end if
    call cpu_time(finish)

    if (debug_mode) then
        call print_matrix_for_method(matprod, dims(1), dims(2), mat_mul_method)
    end if
    print "('Elapsed time = ', es16.10)", finish - start

    ! free memory
    deallocate(matrixA, matrixB, matprod)

end program
