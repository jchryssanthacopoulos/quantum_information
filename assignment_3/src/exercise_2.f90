! Exercise 2
! ==========
!


program exercise_2
    integer n
    complex*16, dimension(:, :), allocatable :: A

    ! matrix size
    n = 5

    ! allocate memory
    allocate(A(n, n))

    ! generatre random hermitian matrix
    A = rand_hermitian_matrix(n)

    ! display matrix
    call print_matrix(A, n)

    deallocate(A)

contains
    function rand_hermitian_matrix(n) result(H)
        integer n
        integer ii, jj
        complex*16, dimension(n, n) :: H

        do ii = 1, n
            do jj = 1, ii
                if (ii /= jj) then
                    ! sample random complex numbers off the diagonal
                    h(ii, jj) = cmplx(RAND(0)*2 - 1, RAND(0)*2 - 1)
                    h(jj, ii) = conjg(h(ii, jj))
                else
                    ! sample real numbers on the diagonal
                    h(ii, ii) = RAND(0)*2 - 1
                end if
            end do
        end do
    end function

    subroutine print_matrix(M, n)
        implicit none

        integer n
        integer ii
        complex*16 M(n, n)

        do ii = 1, n
            print '(*(sp, f7.4, 1x, f7.4, "i", 3x))', M(ii, :)
        end do
    end subroutine

end program
