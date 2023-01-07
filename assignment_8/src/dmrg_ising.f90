! Density Matrix Renormalization Group on Ising Model
! ===================================================
!
! This program computes the ground state energy of the 1D Ising model using the infinte density matrix renormalization 
!   group algorithm
!
! Command-line arguments:
!   N (integer): Number of spin sites
!   lambda (real): Coupling between neighboring sites
!   max_iter (integer): Maximum number of iterations to run
!   thres (real): Threshold for convergence in terms of successive differences in normalized ground state energy
!   debug (logical): Whether to display debug information
!
! Returns:
!   Print the ground state energy, number of iterations run, and whether algorithm converged to standard out
!


module density_matrix_rg
    use ising_hamiltonian
    use density_matrix
    use lin_alg_utils

contains
    ! Run an iteration of density matrix renormalization group starting from subsystem and interaction Hamiltonians
    !
    ! Inputs:
    !   N (integer): Number of sites in each bipartition
    !   H_1 (complex*16 matrix): Hamiltonian of each bipartition
    !   A (complex*16 matrix): Part of interaction Hamiltonian corresponding to left bipartition
    !   lambda (real*8): Interaction with external magnetic field
    !   eigenvalues (real*8 array): Array to store energy eigenvalues
    !   debug (logical): Whether to display debug information
    !
    ! Returns:
    !   Updates H_1 and A in-place
    !
    subroutine iter_density_matrix_rg(N, H_1, A, lambda, eigenvalues, debug)
        implicit none

        integer N, m, d, ndim, ii
        logical debug

        real*8 lambda

        ! Pauli matrices
        complex*16 sigma_x(2, 2), sigma_z(2, 2)

        complex*16, dimension(:, :) :: H_1, A
        complex*16, dimension(:, :), allocatable :: H_enlarged_1, H_12, H_23
        complex*16, dimension(:, :), allocatable :: H

        ! eigenvalues and eigenvectors of H and rho
        real*8, dimension(:), allocatable :: eigenvalues, eigenvalues_rho
        complex*16, dimension(:, :), allocatable :: eigenvectors, eigenvectors_rho

        ! rho matrices
        complex*16, dimension(:, :), allocatable :: rho, rho_reduced_L

        ! projection operators
        complex*16, dimension(:, :), allocatable :: P, P_transpose

        d = 2
        m = d ** N
        ndim = m * d

        sigma_x = get_sigma_x()
        sigma_z = get_sigma_z()

        allocate(H(ndim ** 2, ndim ** 2))
        allocate(H_enlarged_1(ndim, ndim))
        allocate(H_12(ndim, ndim))
        allocate(H_23(ndim ** 2, ndim ** 2))

        allocate(eigenvalues_rho(ndim))
        allocate(eigenvectors(ndim ** 2, ndim ** 2), eigenvectors_rho(ndim, ndim))

        allocate(rho(ndim ** 2, ndim ** 2), rho_reduced_L(ndim, ndim))
        allocate(P(ndim, m), P_transpose(m, ndim))

        ! enlarge the block
        H_12 = tensor_product(A, sigma_x)
        H_enlarged_1 = tensor_product(H_1, identity(1)) + lambda * tensor_product(identity(N), sigma_z) + H_12

        ! interaction (always the same?)
        H_23 = tensor_product(tensor_product(identity(N), sigma_x), tensor_product(sigma_x, identity(N)))

        ! Hamiltonian of left and right enlarged blocks
        H = tensor_product(H_enlarged_1, identity(N + 1)) + tensor_product(identity(N + 1), H_enlarged_1) + H_23

        if (debug) then
            print *, "H_1 = "
            call print_complex_matrix(H_1)
            print *, "A = "
            call print_complex_matrix(A)
            print *, "H_12 = "
            call print_complex_matrix(H_12)
            print *, "H_23 = "
            call print_complex_matrix(H_23)
            print *, "H = "
            call print_complex_matrix(H)
        end if

        ! diagonalize H
        call find_eigenvalues_using_zheev(H, eigenvalues, eigenvectors)

        if (debug) then
            print *, "eigenvalues = ", eigenvalues
        end if

        ! compute density matrix of ground state
        rho = compute_density_matrix(eigenvectors(:, 1))

        ! compute left reduced density matrix
        rho_reduced_L = compute_left_reduced_density_matrix(rho, d, 2 * N + 2, N + 1)

        if (debug) then
            print *, "rho = "
            call print_complex_matrix(rho)
            print *, "rho_reduced_L = "
            call print_complex_matrix(rho_reduced_L)
        end if

        ! diagonalize reduced density matrix
        call find_eigenvalues_using_zheev(rho_reduced_L, eigenvalues_rho, eigenvectors_rho)

        if (debug) then
            print *, "eigenvalues of rho reduced L = ", eigenvalues_rho
        end if

        ! truncate to get projection operator
        do ii = 1, m
            if (debug) then
                print *, "Using eigenvector with eigenvalue = ", eigenvalues_rho(ndim + 1 - ii)
            end if
            P(:, ii) = eigenvectors_rho(:, ndim + 1 - ii)
        end do

        if (debug) then
            print *, "P = "
            call print_complex_matrix(P)
        end if

        P_transpose = transpose(conjg(P))

        ! project into lower-dimensional subspace
        H_1 = matmul(matmul(P_transpose, H_enlarged_1), P)
        A = matmul(matmul(P_transpose, tensor_product(A, identity(1))), P)

        if (debug) then
            print *, "H_1 after iteration = "
            call print_complex_matrix(H_1)
            print *, "A after iteration = "
            call print_complex_matrix(A)
        end if

        deallocate(H, H_enlarged_1, H_12, H_23)
        deallocate(eigenvalues_rho, eigenvectors, eigenvectors_rho, rho, rho_reduced_L, P, P_transpose)

    end subroutine

end module


program dmrg_ising
    use arg_parse
    use density_matrix_rg
    implicit none

    integer d, iter

    ! Pauli matrices
    complex*16 sigma_x(2, 2)

    ! ground-state energy
    real*8 gs_energy, prev_gs_energy

    ! matrices needed to compute Hamiltonian
    complex*16, dimension(:, :), allocatable :: H_1, A

    ! to store energy eigenvalues and eigenvectors
    real*8, dimension(:), allocatable :: eigenvalues, eigenvalues_H
    complex*16, dimension(:, :), allocatable :: eigenvectors_H

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
    write (arg_char, "(l1)") debug
    print *, "debug = ", adjustl(arg_char)

    d = 2

    allocate(eigenvalues(d ** (2 * N + 2)))
    allocate(eigenvalues_H(d ** N))
    allocate(eigenvectors_H(d ** N, d ** N))

    ! compute initial Hamiltonian
    H_1 = lambda * non_interacting_hamiltonian(N) + interacting_hamiltonian(N)

    ! compute operators in interaction Hamiltonian
    sigma_x = get_sigma_x()
    A = tensor_product(identity(N - 1), sigma_x)

    iter = 1
    gs_energy = -1
    prev_gs_energy = 0

    do while ((iter .le. max_iter) .and. abs(gs_energy - prev_gs_energy) > thres)
        prev_gs_energy = gs_energy

        call iter_density_matrix_rg(N, H_1, A, lambda, eigenvalues, debug)

        ! compute energy density
        gs_energy = eigenvalues(1) / dble(d * (N - 1 + iter))

        if (debug) then
            print "('Iteration = ', i3, ', ground state energy = ', f9.4)", iter, gs_energy
        end if

        iter = iter + 1
    end do

    print "('Ground state energy = ', f9.4)", gs_energy
    print "('Max iterations run = ', i4)", iter - 1
    print "('Did converge = ', l1)", abs(gs_energy - prev_gs_energy) <= thres

    deallocate(H_1, A, eigenvalues, eigenvalues_H, eigenvectors_H)

end program
