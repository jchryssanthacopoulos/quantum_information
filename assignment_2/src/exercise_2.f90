! Exercise 2
! ==========
!
! This program implements two different ways of performing matrix multiplication and compares their results and
!   performance against the builtin method, matmul
!
! Example usage:
!  $ exercise_2
! Enter number of rows, columns, and inner dimension:
!  100 100 100
!  Elapsed time for matmul = 1.1700000000E-04
!  Max abs error for row-col-inner = 3.1974423109E-14
!  Elapsed time for row-col-inner = 3.6200000000E-03
!  Max abs error for inner-col-row = 3.1974423109E-14
!  Elapsed time for inner-col-row = 5.0080000000E-03
!
! If you provide non-integer inputs, you'll get the following error:
!  $ exercise_2
! Enter number of rows, columns, and inner dimension:
!  a b c
!  Dimensions need to be integers!
!

#define INT  0
#define REAL 1

#define DTYPE INT

#ifndef DTYPE
#   define DTYPE REAL
#endif

! maximum integer to sample when generating integer matrix
#define MAX_INT_RANGE 100


module exercise_2_utils
    implicit none

    character(len=20) char_input(3)
    integer*4 dims(3)
    integer status
    logical verbose

contains
    ! parse command-line arguments
    subroutine parse_cmd_args()
        implicit none

        integer ii
        character(len=32) arg

        do ii = 1, command_argument_count()
            call get_command_argument(ii, arg)
            select case (arg)
                case ('-v', '--verbose')
                    verbose = .true.
            end select
        end do
    end subroutine

    ! initialize matrix with random data depending on data type
    subroutine init_matrix(M, nrows, ncols)
        implicit none

        integer*4 nrows, ncols

#if DTYPE == REAL
        real*8, intent(inout) :: M(nrows, ncols)
        call random_number(M)
#elif DTYPE == INT
        integer*4, intent(inout) :: M(nrows, ncols)
        integer*4 ii, jj
        real*8 rand_real

        do ii = 1, nrows
            do jj = 1, ncols
                call random_number(rand_real)
                M(ii, jj) = floor((MAX_INT_RANGE + 1)*rand_real)
            end do
        end do
#endif

    end subroutine

    ! print matrix in a nice format
    subroutine print_matrix(M, nrows, ncols)
        implicit none

        integer*4 ii, nrows, ncols

#if DTYPE == REAL
        real*8 M(nrows, ncols)
        do ii = 1, nrows
            print '(20f7.2)', M(ii, 1:ncols)
        end do
#elif DTYPE == INT
        integer*4 M(nrows, ncols)
        do ii = 1, nrows
            print '(20i13)', M(ii, 1:ncols)
        end do
#endif

    end subroutine

    ! display maximum absolute error between two matrices
    function max_abs_error(A, B, nrows, ncols)
        integer*4 nrows, ncols, ii, jj, cnt
        real*8 error, max_abs_error

#if DTYPE == REAL
        real*8 A(nrows, ncols), B(nrows, ncols)
#elif DTYPE == INT
        integer*4 A(nrows, ncols), B(nrows, ncols)
#endif

        cnt = 0

        ! iterate over all entries, computing the absolute difference
        do ii = 1, nrows
            do jj = 1, ncols
                error = abs(A(ii, jj) - B(ii, jj))
                if (cnt == 0) then
                    max_abs_error = error
                    cnt = cnt + 1
                else if (error > max_abs_error) then
                    max_abs_error = error
                end if
            end do
        end do
    end function

    ! check that the input arguments are integers that can be used as matrix dimensions
    ! if conversion fails for any input argument, a non-zero status is returned
    subroutine check_input_args(str, int, status)
        implicit none

        integer*4 ii
        character(len=20) :: str(3)
        integer, intent(out) :: int(3)
        integer, intent(out) :: status

        do ii = 1, 3
            call str2int(str(ii), int(ii), status)
            if (status /= 0) then
                return
            end if
        end do
    end subroutine

    ! converts a string into an integer, returning a non-zero status code if conversion failed
    ! adapted from code found here:
    ! https://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90
    subroutine str2int(str, int, status)
        implicit none

        character(len=*), intent(in) :: str
        integer, intent(out) :: int
        integer, intent(out) :: status

        ! read string into an integer, saving the resulting status
        read(str, *, iostat=status) int
    end subroutine

end module


program exercise_2
    use exercise_2_utils
    implicit none

    integer*4 ii, jj, kk
    real*8 start, finish

#if DTYPE == REAL
    real*8, dimension(:, :), allocatable :: matrixA, matrixB, matprod1, matprod2, matprod3
#elif DTYPE == INT
    integer*4, dimension(:, :), allocatable :: matrixA, matrixB, matprod1, matprod2, matprod3
#endif

    ! read matrix input sizes
    print *, 'Enter number of rows, columns, and inner dimension:'
    read (*, *) char_input

    ! error checking
    call check_input_args(char_input, dims, status)
    if (status /= 0) then
        print *, "Dimensions need to be integers!"
        return
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
    call init_matrix(matrixA, dims(1), dims(3))
    call init_matrix(matrixB, dims(3), dims(2))

    if (verbose) then
        print *, "Running in verbose mode ..."
        print *, "Matrix A = "
        call print_matrix(matrixA, dims(1), dims(3))
        print *, "Matrix B = "
        call print_matrix(matrixB, dims(3), dims(2))
    end if

    ! call intrinsic method
    call cpu_time(start)
        matprod1 = matmul(matrixA, matrixB)
    call cpu_time(finish)

    if (verbose) then
        print *, "Product using matmul = "
        call print_matrix(matprod1, dims(1), dims(2))
    end if
    print "('Elapsed time for matmul = ', es16.10)", finish - start

    ! initialize output matrices
    matprod2 = 0.0
    matprod3 = 0.0

    ! multiply matrices using two different loop orders

    call cpu_time(start)
        do ii = 1, dims(1)
            do jj = 1, dims(2)
                do kk = 1, dims(3)
                    matprod2(ii, jj) = matprod2(ii, jj) + matrixA(ii, kk) * matrixB(kk, jj)
                end do
            end do
        end do
    call cpu_time(finish)

    if (verbose) then
        print *, "Matrix using row-col-inner = "
        call print_matrix(matprod2, dims(1), dims(2))
    end if
    print "('Max abs error for row-col-inner = ', es16.10)", max_abs_error(matprod1, matprod2, dims(1), dims(2))
    print "('Elapsed time for row-col-inner = ', es16.10)", finish - start

    call cpu_time(start)
        do kk = 1, dims(3)
            do jj = 1, dims(2)
                do ii = 1, dims(1)
                    matprod3(ii, jj) = matprod3(ii, jj) + matrixA(ii, kk) * matrixB(kk, jj)
                end do
            end do
        end do
    call cpu_time(finish)

    if (verbose) then
        print *, "Matrix using inner-col-row = "
        call print_matrix(matprod3, dims(1), dims(2))
    end if
    print "('Max abs error for inner-col-row = ', es16.10)", max_abs_error(matprod1, matprod3, dims(1), dims(2))
    print "('Elapsed time for inner-col-row = ', es16.10)", finish - start

    ! free memory
    deallocate(matrixA, matrixB, matprod1, matprod2, matprod3)

end program
