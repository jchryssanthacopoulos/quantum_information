! this module contains utilities related to linear algebra like routines for diagonalization
module lin_alg_utils

contains
    ! Compute first num_eig number of eigenvalues and eigenvectors
    !
    ! Inputs:
    !   A (complex*16 matrix): Matrix to diagonalize
    !   num_eig (integer): Number of eigenvalues to compute
    !   eigenvalues (real*8 array): Array to store eigenvalues
    !   eigenvectors (real*8 matrix): Matrix to store eigenvectors
    !
    subroutine find_k_eigenvalues(A, num_eig, eigenvalues, eigenvectors)
        implicit none

        complex*16, dimension(:, :), allocatable, intent(in) :: A
        real*8, dimension(:), allocatable, intent(inout) :: eigenvalues
        real*8, dimension(:, :), allocatable, intent(inout) :: eigenvectors

        integer N
        integer num_eig

        integer lwork, liwork, info, M, lwmax
        real*8, dimension(:), allocatable :: work
        integer, dimension(:), allocatable :: isuppz, iwork
        real*8, dimension(:, :), allocatable :: supp

        N = size(A, 1)

        lwmax = 100000

        allocate(work(lwmax))
        allocate(iwork(lwmax))
        allocate(isuppz(2 * num_eig))
        allocate(supp(N, N))

        supp = real(A)

        ! compute optimal size of workspace
        lwork = -1
        liwork = -1
        call dsyevr('V', 'I', 'U', N, supp, N, 0.0, 0.0, 1, num_eig, &
            0.0, M, eigenvalues, eigenvectors, N, isuppz, &
            work, lwork, iwork, liwork, info)

        lwork = min(lwmax, int(work(1)))
        liwork = min(lwmax, int(iwork(1)))

        call dsyevr('V', 'I', 'U', N, supp, N, 0.0, 0.0, 1, num_eig, &
            0.0, M, eigenvalues, eigenvectors, N, isuppz, &
            work, lwork, iwork, liwork, info)

        if (info .ne. 0) then
            print *, "Failed to diagonalize"
            stop
        end if

        deallocate(work, iwork, isuppz, supp)

    end subroutine

    ! Compute all eigenvalues and eigenvectors of matrix
    !
    ! Inputs:
    !   A (complex*16 matrix): Matrix to diagonalize
    !   eigenvalues (real*8 array): Array to store eigenvalues
    !   eigenvectors (real*8 matrix): Matrix to store eigenvectors
    !
    subroutine find_eigenvalues_using_zheev(A, eigenvalues, eigenvectors)
        implicit none

        integer ndim, lwork, info
        complex*16, dimension(:, :) :: A
        real*8, dimension(:) :: eigenvalues
        complex*16, dimension(:, :) :: eigenvectors

        complex*16, dimension(:), allocatable :: work
        real*8, dimension(:), allocatable :: rwork

        ndim = size(A, 1)

        lwork = max(1, 2 * ndim - 1)

        allocate(work(max(1, lwork)))
        allocate(rwork(max(1, 3 * ndim - 2)))

        eigenvectors = A

        call zheev("V", "U", ndim, eigenvectors, ndim, eigenvalues, work, lwork, rwork, info)

        if (info .ne. 0) then
            print *, "Failed to diagonalize"
            stop
        end if

        deallocate(work, rwork)

    end subroutine

end module
