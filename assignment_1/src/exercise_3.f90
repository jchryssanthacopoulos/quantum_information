program exercise_3
    implicit none

    integer*4 n_1, n_2, n_3
    integer*4 i, j, k
    real*8 start, finish
    real*8, dimension(:, :), allocatable :: matrix1, matrix2, matrix3, matrix4, matrix5

    ! read matrix input sizes
    print *, 'Enter number of rows, columns, and inner dimension:'
    read (*, *) n_1, n_2, n_3

    ! allocate memory for matrices
    allocate(matrix1(n_1, n_3))
    allocate(matrix2(n_3, n_2))
    allocate(matrix3(n_1, n_2))
    allocate(matrix4(n_1, n_2))
    allocate(matrix5(n_1, n_2))

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
        do i = 1, n_1
            do j = 1, n_2
                do k = 1, n_3
                    matrix3(i, j) = matrix3(i, j) + matrix1(i, k) * matrix2(k, j)
                end do
            end do
        end do
    call cpu_time(finish)
    print "('Max abs error for row-col-inner = ', es16.10)", max_abs_error(matrix3, matrix5, n_1, n_2)
    print "('Elapsed time for row-col-inner = ', es16.10)", finish - start

    call cpu_time(start)
        do k = 1, n_3
            do j = 1, n_2
                do i = 1, n_1
                    matrix4(i, j) = matrix4(i, j) + matrix1(i, k) * matrix2(k, j)
                end do
            end do
        end do
    call cpu_time(finish)
    print "('Max abs error for inner-col-row = ', es16.10)", max_abs_error(matrix4, matrix5, n_1, n_2)
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

end program
