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

        sigma_z = (0d0, 0d0)
        sigma_z(1, 1) = (1d0, 0d0)
        sigma_z(2, 2) = (-1d0, 0d0)

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

        sigma_x = (0d0, 0d0)
        sigma_x(1, 2) = (1d0, 0d0)
        sigma_x(2, 1) = (1d0, 0d0)

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

end module


program ising_model
    use arg_parse
    use ising_hamiltonian
    implicit none

    integer ndim

    ! Hamiltonian
    complex*16, dimension(:, :), allocatable :: H

    ! for computing eigenvalues
    integer lwork, info
    real*8, dimension(:), allocatable :: eigvals
    complex*16, dimension(:), allocatable :: work
    real*8, dimension(:), allocatable :: rwork

    ! variables to clock time to diagonalize
    real*8 start, finish

    ! read number of subsystems and coupling constant
    call parse_cmd_args()
    write (arg_char, "(i8)") N
    print *, "N = ", adjustl(arg_char)
    write (arg_char, "(f7.3)") lambda
    print *, "lambda = ", adjustl(arg_char)
    write (arg_char, "(l1)") debug
    print *, "debug = ", adjustl(arg_char)
    print *, "output_filename = ", output_filename

    ndim = 2 ** N

    lwork = max(1, 2 * ndim - 1)

    allocate(eigvals(ndim))
    allocate(work(max(1, lwork)))
    allocate(rwork(max(1, 3 * ndim - 2)))

    ! compute Hamiltonian
    call cpu_time(start)
        H = lambda * non_interacting_hamiltonian(N)
        H = H + interacting_hamiltonian(N)
    call cpu_time(finish)
    print "('Elapsed time to generate Hamiltonian = ', es16.10)", finish - start

    if (debug) then
        print *, "H is "
        call print_complex_matrix(H)
    end if

    ! compute eigenvalues
    call cpu_time(start)
        call zheev("N", "U", ndim, H, ndim, eigvals, work, lwork, rwork, info)
    call cpu_time(finish)
    print "('Elapsed time to diagonalize = ', es16.10)", finish - start

    open(1, file=output_filename)
    write(1, *) "Energy eigenvalues = "
    write(1, *) eigvals

    deallocate(H, eigvals, work, rwork)

end program
