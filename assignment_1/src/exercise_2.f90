module different_sums
contains
    subroutine sum_int_2(x, y)          
        integer*2 x, y
        integer*2 sum
        sum = x + y
        print *, 'The sum of ', x, 'and', y, 'using INTEGER*2 is', sum
    end subroutine sum_int_2

    subroutine sum_int_4(x, y)          
        integer*4 x, y
        integer*4 sum
        sum = x + y
        print *, 'The sum of ', x, 'and', y, 'using INTEGER*4 is', sum
    end subroutine sum_int_4

end module different_sums


program exercise_2
    use different_sums

    integer n_1, n_2
    integer*2 x, y
    integer*4 x_2, y_2
    real*4 num_1, num_2, sum_1
    real*8 pi

    print *, 'Enter number 1:'
    read (*,*) n_1
    print *, 'Enter number 2:'
    read (*,*) n_2

    ! sum using integer*2
    x = n_1
    y = n_2
    call sum_int_2(x, y)

    ! sum using integer*4
    x_2 = n_1
    y_2 = n_2
    call sum_int_4(x_2, y_2)

    pi = 4.D0*DATAN(1.D0)
    num_1 = pi * 1.0D32
    num_2 = 1.414 * 1D21
    sum_1 = num_1 + num_2
    print *, pi, num_1, num_2, sum_1

end program
