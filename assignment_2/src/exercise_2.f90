! Exercise 2
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


program exercise_2
    use arg_parse
    use debug
    use mat_ops
    implicit none

    character(len=20) char_input(3)
    integer*4 dims(3)
    real*8 start, finish
    real*8, dimension(:, :), allocatable :: matrixA, matrixB, matprod1, matprod2, matprod3

    ! read matrix input sizes
    print *, 'Enter number of rows, columns, and inner dimension:'
    read *, char_input

    ! check if dimensions are integers
    call check_dims_integers(char_input, dims, 3)
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

    ! parse command-line options
    call parse_cmd_args()

    ! allocate memory for matrices
    allocate(matrixA(dims(1), dims(3)))
    allocate(matrixB(dims(3), dims(2)))
    allocate(matprod1(dims(1), dims(2)))
    allocate(matprod2(dims(1), dims(2)))
    allocate(matprod3(dims(1), dims(2)))

    ! initialize input matrices with random data
    call random_number(matrixA)
    call random_number(matrixB)

    if (debug_mode) then
        print *, "Running in debug mode ..."
        call print_input_matrices(matrixA, matrixB, dims)
    end if

    ! call intrinsic method
    call cpu_time(start)
        matprod1 = matmul(matrixA, matrixB)
    call cpu_time(finish)

    if (debug_mode) then
        call print_matrix_for_method(matprod1, dims(1), dims(2), "matmul")
    end if
    print "('Elapsed time for matmul = ', es16.10)", finish - start

    ! multiply matrices using two different loop orders

    matprod2 = 0.0
    call cpu_time(start)
        call matmul_row_col(matrixA, matrixB, matprod2, dims)
    call cpu_time(finish)

    if (debug_mode) then
        call print_matrix_for_method(matprod2, dims(1), dims(2), "row-col")
    end if
    print "('Max abs error for row-col = ', es16.10)", max_abs_error(matprod1, matprod2, dims(1), dims(2))
    print "('Elapsed time for row-col = ', es16.10)", finish - start

    matprod3 = 0.0
    call cpu_time(start)
        call matmul_col_row(matrixA, matrixB, matprod3, dims)
    call cpu_time(finish)

    if (debug_mode) then
        call print_matrix_for_method(matprod3, dims(1), dims(2), "col-row")
    end if
    print "('Max abs error for col-row = ', es16.10)", max_abs_error(matprod1, matprod3, dims(1), dims(2))
    print "('Elapsed time for col-row = ', es16.10)", finish - start

    ! free memory
    deallocate(matrixA, matrixB, matprod1, matprod2, matprod3)

end program
