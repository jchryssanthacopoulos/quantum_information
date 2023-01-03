! this module contains functions for computing the Hamiltonian of the 1D Ising model
module ising_hamiltonian

contains
    ! compute the tensor product between two complex matrices
    !
    ! Inputs:
    !   A (complex*16 matrix): First matrix
    !   B (complex*16 matrix): Second matrix
    !
    ! Returns:
    !   TP (complex*16 matrix): Tensor product
    !
    function tensor_product(A, B) result(TP)
        implicit none

        integer ii, jj, kk, mm
        integer idx1, idx2
        integer N_A(2), N_B(2), N_TP(2)

        complex*16, dimension(:, :) :: A, B
        complex*16, dimension(:, :), allocatable :: TP

        N_A = shape(A)
        N_B = shape(B)
        N_TP = (/N_A(1) * N_B(1), N_A(2) * N_B(2)/)

        allocate(TP(N_TP(1), N_TP(2)))

        do ii = 1, N_A(1)
            do jj = 1, N_A(2)
                do kk = 1, N_B(1)
                    do mm = 1, N_B(2)
                        idx1 = (ii - 1) * N_B(1) + kk
                        idx2 = (jj - 1) * N_B(2) + mm
                        TP(idx1, idx2) = A(ii, jj) * B(kk, mm)
                    end do
                end do
            end do
        end do

    end function

    ! compute the identity matrix for an N-qubit system
    !
    ! Inputs:
    !   N (integer): Number of qubits
    !
    ! Returns:
    !   I (complex*16 matrix): Identity matrix with dimension 2 ** N
    !
    function identity(N) result(I)
        implicit none

        integer N
        integer dim
        integer ii
        complex*16, dimension(:, :), allocatable :: I

        dim = 2 ** N

        allocate(I(dim, dim))

        I = (0d0, 0d0)

        do ii = 1, dim
            I(ii, ii) = (1d0, 0d0)
        end do

    end function

    ! compute non-interacting part of Hamiltonian
    !
    ! Inputs:
    !   N (integer): Number of spin sites
    !
    ! Returns:
    !   H_0 (complex*16 matrix): Non-interacting Hamiltonian
    !
    function non_interacting_hamiltonian(N) result(H_0)
        implicit none

        integer N, dim
        integer ii
        complex*16 sigma_z(2, 2)
        complex*16, dimension(:, :), allocatable :: H_0, H_0_i

        complex*16, dimension(:, :), allocatable :: I1, I2, prod1

        sigma_z = get_sigma_z()

        dim = size(sigma_z, 1) ** N

        allocate(H_0(dim, dim))

        H_0 = (0d0, 0d0)

        do ii = 1, N
            I1 = identity(ii - 1)
            I2 = identity(N - ii)
            prod1 = tensor_product(I1, sigma_z)
            H_0_i = tensor_product(prod1, I2)

            H_0 = H_0 + H_0_i

            deallocate(I1, I2, prod1, H_0_i)
        end do

    end function

    ! compute interacting part of Hamiltonian
    !
    ! Inputs:
    !   N (integer): Number of spin sites
    !
    ! Returns:
    !   H_int (complex*16 matrix): Interacting Hamiltonian
    !
    function interacting_hamiltonian(N) result(H_int)
        implicit none

        integer N, dim
        integer ii
        complex*16 sigma_x(2, 2)
        complex*16, dimension(:, :), allocatable :: H_int, H_int_i

        complex*16, dimension(:, :), allocatable :: I1, I2, prod1, prod2

        sigma_x = get_sigma_x()

        dim = size(sigma_x, 1) ** N

        allocate(H_int(dim, dim))

        H_int = (0d0, 0d0)

        do ii = 1, N - 1
            I1 = identity(ii - 1)
            I2 = identity(N - ii - 1)
            prod1 = tensor_product(I1, sigma_x)
            prod2 = tensor_product(prod1, sigma_x)
            H_int_i = tensor_product(prod2, I2)

            H_int = H_int + H_int_i

            deallocate(I1, I2, prod1, prod2, H_int_i)
        end do

    end function

    ! return sigma_x Pauli matrix
    !
    ! Returns:
    !   sigma_x (complex*16 matrix): x Pauli matrix
    !
    function get_sigma_x() result(sigma_x)
        complex*16 sigma_x(2, 2)

        sigma_x = (0d0, 0d0)
        sigma_x(1, 2) = (1d0, 0d0)
        sigma_x(2, 1) = (1d0, 0d0)
    end function

    ! return sigma_z Pauli matrix
    !
    ! Returns:
    !   sigma_z (complex*16 matrix): z Pauli matrix
    !
    function get_sigma_z() result(sigma_z)
        complex*16 sigma_z(2, 2)

        sigma_z = (0d0, 0d0)
        sigma_z(1, 1) = (1d0, 0d0)
        sigma_z(2, 2) = (-1d0, 0d0)
    end function

    ! print square complex matrix in a nice format
    !
    ! Inputs:
    !   M (complex*16 matrix): Matrix to print
    !
    subroutine print_complex_matrix(M)
        implicit none

        integer ii
        complex*16, dimension(:, :) :: M

        do ii = 1, size(M, 1)
            print '(*(sp, f8.5, 1x, f8.5, "i", 3x))', M(ii, :)
        end do
    end subroutine

    ! print matrix in a nice format
    !
    ! Inputs:
    !   M (real*8 matrix): Matrix to print
    !   nrows (integer): Number of rows
    !   ncols (integer): Number of columns
    !
    subroutine print_matrix(M, nrows, ncols)
        implicit none

        integer*4 ii, nrows, ncols
        real*8 M(nrows, ncols)

        do ii = 1, nrows
            print '(20f7.2)', M(ii, 1:ncols)
        end do
    end subroutine

end module
