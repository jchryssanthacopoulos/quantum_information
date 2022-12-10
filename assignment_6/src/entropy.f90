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
    function compute_entropy(rho) result(S)
        implicit none

        integer ii

        ! density matrix
        integer ndim
        complex*16, dimension(:, :) :: rho

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

        allocate(eigvals(ndim))
        allocate(work(max(1, lwork)))
        allocate(rwork(max(1, 3 * ndim - 2)))

        ! compute eigenvalues
        call zheev("N", "U", ndim, rho, ndim, eigvals, work, lwork, rwork, info)

        print *, "Eigenvalues = ", eigvals

        S = 0
        do ii = 1, ndim
            if (abs(eigvals(ii)) .gt. nonzero_eigval_threshold) then
                S = S + eigvals(ii) * log(eigvals(ii))
            end if
        end do

        deallocate(eigvals, work, rwork)

    end function

end module
