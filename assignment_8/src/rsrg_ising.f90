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
!   output_filename (character): Name of file to save solution
!   debug (integer): Whether to display debug information
!
! Returns:
!   Saves the energy and associated metadata to the specified file
!


module real_space_rg
    use ising_hamiltonian

contains
    ! Run an iteration of real-space renormalization group starting from subsystem and interaction Hamiltonians
    !
    ! Inputs:
    !   N (integer): Number of sites in each bipartition
    !   H_N (complex*16 matrix): Hamiltonian of each bipartition
    !   A (complex*16 matrix): Part of interaction Hamiltonian corresponding to left bipartition
    !   B (complex*16 matrix): Part of interaction Hamiltonian corresponding to right bipartition
    !   eigenvalues (real*8 array): Array to store energy eigenvalues
    !
    ! Returns:
    !   Updates H_N, A, and B in-place
    !
    subroutine iter_real_space_rg(N, H_N, A, B, eigenvalues)
        implicit none

        integer N
        integer ndim

        complex*16, dimension(:, :) :: H_N, A, B
        real*8, dimension(:), allocatable :: eigenvalues

        complex*16, dimension(:, :), allocatable :: H_2N, A_2N, B_2N
        real*8, dimension(:, :), allocatable :: P, P_transpose

        ndim = 2 ** N

        allocate(H_2N(ndim ** 2, ndim ** 2))
        allocate(A_2N(ndim ** 2, ndim ** 2))
        allocate(B_2N(ndim ** 2, ndim ** 2))
        allocate(P(ndim ** 2, ndim))
        allocate(P_transpose(ndim, ndim ** 2))

        H_2N = tensor_product(H_N, identity(N)) + tensor_product(identity(N), H_N) + tensor_product(A, B)

        ! compute first ndim eigenvalues and eigenvectors
        call find_k_eigenvalues(H_2N, ndim, eigenvalues, P)

        P_transpose = transpose(P)

        ! project Hamiltonian into lower-dimensional subspace
        H_N = matmul(matmul(P_transpose, H_2N), P)

        ! project A and B
        A_2N = tensor_product(A, identity(N))
        B_2N = tensor_product(identity(N), B)
        A = matmul(matmul(P_transpose, A_2N), P)
        B = matmul(matmul(P_transpose, B_2N), P)

        deallocate(H_2N, A_2N, B_2N, P, P_transpose)

    end subroutine

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
        LWORK = -1
        LIWORK = -1
        call dsyevr('V', 'I', 'U', N, supp, N, 0.0, 0.0, 1, num_eig, &
            0.0, M, eigenvalues, eigenvectors, N, isuppz, &
            work, lwork, iwork, liwork, info)

        lwork = min(lwmax, int(work(1)))
        liwork = min(lwmax, int(iwork(1)))

        call dsyevr('V', 'I', 'U', N, supp, N, 0.0, 0.0, 1, num_eig, &
            0.0, M, eigenvalues, eigenvectors, N, isuppz, &
            work, lwork, iwork, liwork, info)

        if (info .ne. 0) then
            print *, 'Failed to diagonalize'
            stop
        end if

        deallocate(work, iwork, isuppz, supp)

    end subroutine

end module


program rsrg_ising
    use arg_parse
    use real_space_rg
    implicit none

    integer iter
    integer ndim

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
    write (arg_char, "(l1)") debug
    print *, "debug = ", adjustl(arg_char)
    print *, "output_filename = ", output_filename

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

        call iter_real_space_rg(N, H, A, B, eigenvalues)
        gs_energy = eigenvalues(1) / dble(N)

        print *, "Iteration = ", iter, ", ground-state energy = ", gs_energy

        iter = iter + 1
    end do

    deallocate(H, A, B, eigenvalues)

end program
