program exercise_1_matmul
    integer*4 n_1, n_2, n_3, n_4
    real*8, dimension(:, :), allocatable :: matrix1, matrix2, matrix3, matrix4
    integer*4 :: i, j
    real*8 :: start, finish

    ! read matrix input sizes
    print *, 'Enter matrix dimensions of first matrix:'
    read (*,*) n_1, n_2
    print *, 'Enter matrix dimensions of second matrix:'
    read (*,*) n_3, n_4

    ! return if inner dimensions do not agree
    if (n_2 /= n_3) then
        print *, "Inner dimensions must agree!"
        return
    end if

    ! allocate memory for matrices
    allocate(matrix1(n_1, n_2))
    allocate(matrix2(n_3, n_4))
    allocate(matrix3(n_1, n_4))
    allocate(matrix4(n_1, n_4))

    ! initialize input and output matrices
    call random_number(matrix1)
    call random_number(matrix2)
    matrix3 = 0.0

    call cpu_time(start)
        do i = 1, n_1
            do j = 1, n_4
                do k = 1, n_2
                    matrix3(i, j) = matrix3(i, j) + matrix1(i, k) * matrix2(k, j)
                end do
            end do
        end do
    call cpu_time(finish)
    print '("Time = ",f11.9," seconds.")', finish - start

    call cpu_time(start)
        matrix4 = matmul(matrix1, matrix2)
    call cpu_time(finish)
    print '("Time = ",f11.9," seconds.")', finish - start

    ! display results
    print *, "Matrix A ="
    call print_matrix(matrix1, n_1, n_2)
    print *, "Matrix B ="
    call print_matrix(matrix2, n_3, n_4)
    print *, "Custom: A*B ="
    call print_matrix(matrix3, n_1, n_4)
    print *, "Matmul: A*B ="
    call print_matrix(matrix4, n_1, n_4)

    ! free memory
    deallocate(matrix1, matrix2, matrix3, matrix4)

contains
    subroutine print_matrix(b,n,m)
        integer*4 :: n, m
        real*8 :: b(n, m)
        do i = 1, n
            print '(20f7.2)', b(i,1:m)
        enddo
    endsubroutine

end program exercise_1_matmul
