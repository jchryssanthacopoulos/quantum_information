! Exercise 3
! ==========
!
! This program implements two different ways of performing matrix multiplication and compares their results and
!   and performance against the builtin method, matmul
!
! Example usage:
!   $ exercise_3
! Enter number of rows, columns, and inner dimension:
!  100 100 100
!  Elapsed time for matmul = 1.1700000000E-04
!  Max abs error for row-col-inner = 3.1974423109E-14
!  Elapsed time for row-col-inner = 3.6200000000E-03
!  Max abs error for inner-col-row = 3.1974423109E-14
!  Elapsed time for inner-col-row = 5.0080000000E-03
!
! If you provide non-integer inputs, you'll get the following error:
!   $ exercise_3
! Enter number of rows, columns, and inner dimension:
!  a b c
!  Dimensions need to be integers!
!


program exercise_3
    implicit none

    character(len=20) char_input(3)
    integer*4 dims(3)
    integer status

    integer*4 i, j, k
    real*8 start, finish
    real*8, dimension(:, :), allocatable :: matrix1, matrix2, matrix3, matrix4, matrix5

    ! read matrix input sizes
    print *, 'Enter number of rows, columns, and inner dimension:'
    read (*, *) char_input

    ! error checking
    call check_input_args(char_input, dims, status)
    if (status /= 0) then
        print *, "Dimensions need to be integers!"
        return
    end if

    ! allocate memory for matrices
    allocate(matrix1(dims(1), dims(3)))
    allocate(matrix2(dims(3), dims(2)))
    allocate(matrix3(dims(1), dims(2)))
    allocate(matrix4(dims(1), dims(2)))
    allocate(matrix5(dims(1), dims(2)))

    ! initialize input matrices with random data
    call random_number(matrix1)
    call random_number(matrix2)

    ! call intrinsic method
    call cpu_time(start)
        matrix5 = matmul(matrix1, matrix2)
    call cpu_time(finish)
    print "('Elapsed time for matmul = ', es16.10)", finish - start

    ! initialize output matrices
    matrix3 = 0.0
    matrix4 = 0.0

    ! multiply matrices using two different loop orders

    call cpu_time(start)
        do i = 1, dims(1)
            do j = 1, dims(2)
                do k = 1, dims(3)
                    matrix3(i, j) = matrix3(i, j) + matrix1(i, k) * matrix2(k, j)
                end do
            end do
        end do
    call cpu_time(finish)
    print "('Max abs error for row-col-inner = ', es16.10)", max_abs_error(matrix3, matrix5, dims(1), dims(2))
    print "('Elapsed time for row-col-inner = ', es16.10)", finish - start

    call cpu_time(start)
        do k = 1, dims(3)
            do j = 1, dims(2)
                do i = 1, dims(1)
                    matrix4(i, j) = matrix4(i, j) + matrix1(i, k) * matrix2(k, j)
                end do
            end do
        end do
    call cpu_time(finish)
    print "('Max abs error for inner-col-row = ', es16.10)", max_abs_error(matrix4, matrix5, dims(1), dims(2))
    print "('Elapsed time for inner-col-row = ', es16.10)", finish - start

    ! free memory
    deallocate(matrix1, matrix2, matrix3, matrix4, matrix5)

contains
    ! display maximum absolute error between two matrices
    function max_abs_error(a, b, nrows, ncols)
        integer*4 nrows, ncols, cnt
        real*8 a(nrows, ncols), b(nrows, ncols)
        real*8 error, max_abs_error

        cnt = 0

        ! iterate over all entries, computing the absolute difference
        do i = 1, nrows
            do j = 1, ncols
                error = abs(a(i, j) - b(i, j))
                if (cnt == 0) then
                    max_abs_error = error
                    cnt = cnt + 1
                else if (error > max_abs_error) then
                    max_abs_error = error
                end if
            end do
        end do

        return
    end function

    ! check that the input arguments are integers that can be used as matrix dimensions
    ! if conversion fails for any input argument, a non-zero status is returned
    subroutine check_input_args(str, int, status)
        implicit none

        character(len=20) :: str(3)
        integer, intent(out) :: int(3)
        integer, intent(out) :: status

        do i = 1, 3
            call str2int(str(i), int(i), status)
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

end program
