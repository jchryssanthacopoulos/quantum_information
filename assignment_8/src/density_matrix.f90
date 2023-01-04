! this module contains functions for computing density matrices
module density_matrix

contains
    ! compute density matrix of state
    !
    ! Inputs:
    !   state (complex*16 array): Wavefunction to compute density matrix for
    !
    ! Returns:
    !   rho (complex*16 matrix): Density matrix
    !
    function compute_density_matrix(state) result(rho)
        implicit none

        integer dim
        complex*16, dimension(:) :: state
        complex*16, dimension(:, :), allocatable :: rho
        complex*16, dimension(: ,:), allocatable :: bra, ket

        dim = size(state)

        allocate(rho(dim, dim), bra(1, dim), ket(dim, 1))

        ket(:, 1) = state
        bra(1, :) = conjg(state)

        rho = matmul(ket, bra)

        deallocate(bra, ket)

    end function

    ! compute left reduced density matrix after tracing over M subsystems
    !
    ! Inputs:
    !   rho (complex*16 array): Density matrix with dimensions (D ** N, D ** N)
    !   N (integer): Number of subsystems
    !   D (integer): Dimension of each subsystem
    !   M (integer): Number of subsystems to trace over
    !
    ! Returns:
    !   rho_reduced_L (complex*16 matrix): Left reduced density matrix with dimensions (D ** (N - M), D ** (N - M))
    !
    function compute_left_reduced_density_matrix(rho, D, N, M) result(rho_reduced_L)
        implicit none

        integer D, N, M
        integer dim
        integer debug_level
        integer ii, jj, kk, idx1, idx2

        complex*16, dimension(:, :) :: rho
        complex*16, dimension(:, :), allocatable :: rho_reduced_L

        dim = D ** (N - M)

        allocate(rho_reduced_L(dim, dim))

        rho_reduced_L = 0

        do ii = 1, dim
            do jj = 1, dim
                do kk = 1, D ** M
                    idx1 = kk + (ii - 1) * D ** M
                    idx2 = kk + (jj - 1) * D ** M
                    rho_reduced_L(ii, jj) = rho_reduced_L(ii, jj) + rho(idx1, idx2)
                end do
            end do
        end do

    end function

    ! compute right reduced density matrix after tracing over N - M subsystems
    !
    ! Inputs:
    !   rho (complex*16 array): Density matrix with dimensions (D ** N, D ** N)
    !   N (integer): Number of subsystems
    !   D (integer): Dimension of each subsystem
    !   M (integer): Number of subsystems to get density matrix for
    !
    ! Returns:
    !   rho_reduced_R (complex*16 matrix): Right reduced density matrix with dimensions (D ** M, D ** M)
    !
    function compute_right_reduced_density_matrix(rho, D, N, M) result(rho_reduced_R)
        implicit none

        integer D, N, M
        integer dim
        integer debug_level
        integer ii, jj, kk, idx1, idx2

        complex*16, dimension(:, :) :: rho
        complex*16, dimension(:, :), allocatable :: rho_reduced_R

        dim = D ** M

        allocate(rho_reduced_R(dim, dim))

        rho_reduced_R = 0

        do ii = 1, dim
            do jj = 1, dim
                do kk = 1, D ** (N - M)
                    idx1 = ii + (kk - 1) * D ** M
                    idx2 = jj + (kk - 1) * D ** M
                    rho_reduced_R(ii, jj) = rho_reduced_R(ii, jj) + rho(idx1, idx2)
                end do
            end do
        end do

    end function

end module
