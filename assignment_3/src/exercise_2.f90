! Exercise 2
! ==========
!


program exercise_2
    implicit none

    integer ii
    integer n

    ! to store hermitian matrix
    complex*16, dimension(:, :), allocatable :: A

    ! variables for computing eigenvalues and derived quantities
    integer lwork, info
    real*8 ave_delta_eigvals
    real*8, dimension(:), allocatable :: eigvals, norm_eigval_spacings
    complex*8, dimension(:), allocatable :: work
    real*8, dimension(:), allocatable :: rwork

    ! matrix size
    n = 5

    ! used for computing eigenvalues
    lwork = max(1, 2 * n - 1)

    ! allocate memory
    allocate(A(n, n))
    allocate(eigvals(n))
    allocate(norm_eigval_spacings(n - 1))
    allocate(work(max(1, lwork)))
    allocate(rwork(max(1, 3 * n - 2)))

    ! generatre random hermitian matrix
    A = rand_hermitian_matrix(n)

    ! display matrix
    call print_matrix(A, n)

    ! compute eigenvalues of matrix A (ordered in ASCENDING ORDER)
    call zheev('N', 'U', n, A, n, eigvals, work, lwork, rwork, info)

    ! compute eigenvalue spacings and average
    ave_delta_eigvals = (eigvals(n) - eigvals(1)) / (n - 1)
    do ii = 1, n - 1
        norm_eigval_spacings(ii) = (eigvals(ii + 1) - eigvals(ii)) / ave_delta_eigvals
    end do

    ! print the eigenvalues
    print *, "Info =", info
    print *, "Eigenvalues =", eigvals
    print *, "Normalized eigenvalue spacings =", norm_eigval_spacings
    print *, "Average eigenvalue spacing =", ave_delta_eigvals

    deallocate(A, eigvals, norm_eigval_spacings, work, rwork)

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
