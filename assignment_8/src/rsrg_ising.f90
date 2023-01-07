! Real-space Renormalization Group on Ising Model
! ===============================================
!
! This program computes the ground state energy of the 1D Ising model using the real-space renormalization group
!   algorithm
!
! Command-line arguments:
!   N (integer): Number of spin sites
!   lambda (real): Coupling between neighboring sites
!   max_iter (integer): Maximum number of iterations to run
!   thres (real): Threshold for convergence in terms of successive differences in normalized ground state energy
!   diag_method (character): Method to use to diagonalize (i.e., dsyevr or zheev)
!   debug (logical): Whether to display debug information
!
! Returns:
!   Print the ground state energy, number of iterations run, and whether algorithm converged to standard out
!


module real_space_rg
    use ising_hamiltonian
    use lin_alg_utils

contains
    ! Run an iteration of real-space renormalization group starting from subsystem and interaction Hamiltonians
    !
    ! Inputs:
    !   N (integer): Number of sites in each bipartition
    !   H_N (complex*16 matrix): Hamiltonian of each bipartition
    !   A (complex*16 matrix): Part of interaction Hamiltonian corresponding to left bipartition
    !   B (complex*16 matrix): Part of interaction Hamiltonian corresponding to right bipartition
    !   eigenvalues (real*8 array): Array to store energy eigenvalues
    !   diag_method (character): Method to use to diagonalize (i.e., dsyevr or zheev)
    !   debug (logical): Whether to display debug information
    !
    ! Returns:
    !   Updates H_N, A, and B in-place
    !
    subroutine iter_real_space_rg(N, H_N, A, B, eigenvalues, diag_method, debug)
        implicit none

        integer N, ndim, ii
        logical debug
        character(len=*) diag_method

        complex*16, dimension(:, :) :: H_N, A, B
        complex*16, dimension(:, :), allocatable :: H_2N, A_2N, B_2N
        real*8, dimension(:), allocatable :: eigenvalues

        real*8, dimension(:, :), allocatable :: P, P_transpose
        complex*16, dimension(:, :), allocatable :: full_P, full_P_trunc, full_P_trunc_transpose

        ndim = 2 ** N

        allocate(H_2N(ndim ** 2, ndim ** 2))
        allocate(A_2N(ndim ** 2, ndim ** 2))
        allocate(B_2N(ndim ** 2, ndim ** 2))
        allocate(P(ndim ** 2, ndim))
        allocate(P_transpose(ndim, ndim ** 2))
        allocate(full_P(ndim ** 2, ndim ** 2))
        allocate(full_P_trunc(ndim ** 2, ndim))
        allocate(full_P_trunc_transpose(ndim, ndim ** 2))

        H_2N = tensor_product(H_N, identity(N)) + tensor_product(identity(N), H_N) + tensor_product(A, B)
        A_2N = tensor_product(A, identity(N))
        B_2N = tensor_product(identity(N), B)

        if (debug) then
            print *, "H_N = "
            call print_complex_matrix(H_N)
            print *, "A = "
            call print_complex_matrix(A)
            print *, "B = "
            call print_complex_matrix(B)
            print *, "H_2N = "
            call print_complex_matrix(H_2N)
        end if

        if (diag_method .eq. "dsyevr") then
            ! compute first ndim eigenvalues and eigenvectors
            call find_k_eigenvalues(H_2N, ndim, eigenvalues, P)

            if (debug) then
                print *, "P = "
                call print_matrix(P, ndim ** 2, ndim)
            end if

            P_transpose = transpose(P)

            H_N = matmul(matmul(P_transpose, H_2N), P)
            A = matmul(matmul(P_transpose, A_2N), P)
            B = matmul(matmul(P_transpose, B_2N), P)
        else
            ! compute using zheev
            call find_eigenvalues_using_zheev(H_2N, eigenvalues, full_P)

            ! truncate
            do ii = 1, ndim
                full_P_trunc(:, ii) = full_P(:, ii)
            end do

            if (debug) then
                print *, "P = "
                call print_complex_matrix(full_P_trunc)
            end if

            full_P_trunc_transpose = transpose(conjg(full_P_trunc))

            H_N = matmul(matmul(full_P_trunc_transpose, H_2N), full_P_trunc)
            A = matmul(matmul(full_P_trunc_transpose, A_2N), full_P_trunc)
            B = matmul(matmul(full_P_trunc_transpose, B_2N), full_P_trunc)
        end if

        if (debug) then
            print *, "H_N after iteration = "
            call print_complex_matrix(H_N)
            print *, "A after iteration = "
            call print_complex_matrix(A)
            print *, "B after iteration = "
            call print_complex_matrix(B)
        end if

        deallocate(H_2N, A_2N, B_2N, P, P_transpose, full_P, full_P_trunc, full_P_trunc_transpose)

    end subroutine

end module


program rsrg_ising
    use arg_parse
    use real_space_rg
    implicit none

    integer iter, ndim

    ! Pauli matrices
    complex*16 sigma_x(2, 2)

    ! ground-state energy
    real*8 gs_energy, prev_gs_energy

    ! matrices needed to compute Hamiltonian
    complex*16, dimension(:, :), allocatable :: H, A, B

    ! to store energy eigenvalues
    real*8, dimension(:), allocatable :: eigenvalues

    ! read number of subsystems and coupling constant
    call parse_cmd_args()
    write (arg_char, "(i8)") N
    print *, "N = ", adjustl(arg_char)
    write (arg_char, "(i8)") max_iter
    print *, "max_iter = ", adjustl(arg_char)
    write (arg_char, "(f7.3)") lambda
    print *, "lambda = ", adjustl(arg_char)
    write (arg_char, "(e8.3)") thres
    print *, "thres = ", adjustl(arg_char)
    print *, "diag method = ", diag_method
    write (arg_char, "(l1)") debug
    print *, "debug = ", adjustl(arg_char)

    ndim = 2 ** N

    allocate(eigenvalues(ndim ** 2))

    ! compute initial Hamiltonian
    H = lambda * non_interacting_hamiltonian(N) + interacting_hamiltonian(N)

    ! compute operators in interaction Hamiltonian
    sigma_x = get_sigma_x()
    A = tensor_product(identity(N - 1), sigma_x)
    B = tensor_product(sigma_x, identity(N - 1))

    iter = 1
    gs_energy = -1
    prev_gs_energy = 0

    do while ((iter .le. max_iter) .and. abs(gs_energy - prev_gs_energy) > thres)
        prev_gs_energy = gs_energy

        H = H * 0.5
        A = 1 / sqrt(2.) * A
        B = 1 / sqrt(2.) * B

        call iter_real_space_rg(N, H, A, B, eigenvalues, diag_method, debug)

        ! compute energy density
        gs_energy = eigenvalues(1) / dble(N)

        if (debug) then
            print "('Iteration = ', i3, ', ground state energy = ', f9.4)", iter, gs_energy
        end if

        iter = iter + 1
    end do

    print "('Ground state energy = ', f9.4)", gs_energy
    print "('Max iterations run = ', i3)", iter - 1
    print "('Did converge = ', l1)", abs(gs_energy - prev_gs_energy) <= thres

    deallocate(H, A, B, eigenvalues)

end program
