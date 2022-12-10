! this module computes the entropy associated with a density matrix
module entropy

contains
    ! compute entropy of given density matrix
    !
    ! Inputs:
    !   rho (complex*16 array): Density matrix
    !
    ! Returns:
    !   S (real*8): Entropy
    !
    function compute_entropy(rho, debug) result(S)
        implicit none

        integer ii
        logical debug

        ! density matrix
        integer ndim
        complex*16, dimension(:, :) :: rho
        complex*16, dimension(:, :), allocatable :: rho_copy

        ! for computing eigenvalues
        integer lwork, info
        real*8, dimension(:), allocatable :: eigvals
        complex*16, dimension(:), allocatable :: work
        real*8, dimension(:), allocatable :: rwork

        ! entropy
        real*8, parameter :: nonzero_eigval_threshold = 1E-8
        real*8 S

        ndim = size(rho, 1)

        lwork = max(1, 2 * ndim - 1)

        allocate(rho_copy(ndim, ndim))
        allocate(eigvals(ndim))
        allocate(work(max(1, lwork)))
        allocate(rwork(max(1, 3 * ndim - 2)))

        ! make a copy of the density matrix
        rho_copy = rho

        ! compute eigenvalues
        call zheev("N", "U", ndim, rho_copy, ndim, eigvals, work, lwork, rwork, info)

        if (debug) then
            print *, "Eigenvalues of density matrix = ", eigvals
        end if

        S = 0
        do ii = 1, ndim
            if (abs(eigvals(ii)) .gt. nonzero_eigval_threshold) then
                S = S - eigvals(ii) * log(eigvals(ii))
            end if
        end do

        deallocate(rho_copy, eigvals, work, rwork)

    end function

end module
