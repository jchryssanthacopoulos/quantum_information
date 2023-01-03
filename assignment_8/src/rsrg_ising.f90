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
    subroutine iter_real_space_rg(N, H_N, A, B, eigenvalues, debug)
        implicit none

        integer N
        integer ndim
        integer ii
        logical debug

        complex*16, dimension(:, :) :: H_N, A, B
        real*8, dimension(:), allocatable :: eigenvalues

        complex*16, dimension(:, :), allocatable :: H_2N, A_2N, B_2N
        real*8, dimension(:, :), allocatable :: P, P_transpose

        ! full set of eigenvectors computed using zheev
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

        ! compute first ndim eigenvalues and eigenvectors
        call find_k_eigenvalues(H_2N, ndim, eigenvalues, P)

        ! compute using zheev
        call find_eigenvalues_using_zheev(H_2N, eigenvalues, full_P)

        ! truncate
        do ii = 1, ndim
            full_P_trunc(:, ii) = full_P(:, ii)
        end do

        if (debug) then
            print *, "P = "
            call print_matrix(P, ndim ** 2, ndim)

            print *, "P trunc = "
            call print_complex_matrix(full_P_trunc)
        end if

        P_transpose = transpose(P)
        full_P_trunc_transpose = transpose(conjg(full_P_trunc))

        ! project Hamiltonian into lower-dimensional subspace
        H_N = matmul(matmul(full_P_trunc_transpose, H_2N), full_P_trunc)
        ! H_N = matmul(matmul(P_transpose, H_2N), P)

        ! project A and B
        A_2N = tensor_product(A, identity(N))
        B_2N = tensor_product(identity(N), B)
        A = matmul(matmul(full_P_trunc_transpose, A_2N), full_P_trunc)
        B = matmul(matmul(full_P_trunc_transpose, B_2N), full_P_trunc)
        ! A = matmul(matmul(P_transpose, A_2N), P)
        ! B = matmul(matmul(P_transpose, B_2N), P)

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
            print *, 'Failed to diagonalize'
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

        deallocate(work, rwork)

    end subroutine

    ! Print matrix in a nice format
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

        H = H * 0.5
        A = 1 / sqrt(2.) * A
        B = 1 / sqrt(2.) * B

        call iter_real_space_rg(N, H, A, B, eigenvalues, debug)

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
