! Solve 1D Ising Model
! ====================
!
! This program generates the Hamiltonian of the 1D Ising model and solves for the eigenvalues and eigenvectors
!
! Command-line arguments:
!   N (integer): Number of spin sites
!   lambda (real): Coupling between neighboring sites
!   output_filename (character): Name of file to save solution
!   debug (integer): Whether to display debug information
!
! Returns:
!   Saves the Hamiltonian and solutions to specified file
!


module ising_hamiltonian

contains
    ! Compute the tensor product between two complex matrices
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

    ! Compute the identity matrix for an N-qubit system
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

        do ii = 1, dim
            I(ii, ii) = (1d0, 0d0)
        end do

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
            print '(*(sp, f7.4, 1x, f7.4, "i", 3x))', M(ii, :)
        end do
    end subroutine

end module


program ising_model
    use arg_parse
    use ising_hamiltonian
    implicit none

    complex*16, dimension(:, :), allocatable :: I
    complex*16, dimension(:, :), allocatable :: TP

    ! Pauli matrices
    complex*16 sigma_x(2, 2), sigma_z(2, 2)

    ! read number of subsystems and coupling constant
    call parse_cmd_args()
    write (arg_char, "(i8)") N
    print *, "N = ", adjustl(arg_char)
    write (arg_char, "(f7.3)") lambda
    print *, "lambda = ", adjustl(arg_char)
    write (arg_char, "(l1)") debug
    print *, "debug = ", adjustl(arg_char)
    print *, "output_filename = ", output_filename

    ! define Pauli matrices
    sigma_x = (0d0, 0d0)
    sigma_x(1, 2) = (1d0, 0d0)
    sigma_x(2, 1) = (1d0, 0d0)

    sigma_z = (0d0, 0d0)
    sigma_z(1, 1) = (1d0, 0d0)
    sigma_z(2, 2) = (-1d0, 0d0)

    ! compute tensor product
    I = identity(1)
    TP = tensor_product(I, sigma_x)
    print *, "Tensor product between I and sigma_x is "
    call print_complex_matrix(TP)

    TP = tensor_product(I, sigma_z)
    print *, "Tensor product between I and sigma_z is "
    call print_complex_matrix(TP)

    deallocate(I, TP)

end program
