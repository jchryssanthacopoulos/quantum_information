! Generate density matrix
! =======================
!
! This program generates the density matrix for a many-body quantum state and its left and right reduced matrices
!
! Command-line arguments:
!   N (integer): Number of subsystems in the quantum state
!   dim (integer): Number of dimensions of each subsystem
!   type (character): Type of system to compute (options are "separable", "bell", and "generic")
!   M (integer): Number of subsystems in the right system
!   output_filename (character): Name of file to save density matrices
!   debug (integer): Debug level (options are 0, 1, and 2)
!
! Returns:
!   Saves the density matrices to the specified file
!


module many_body_quantum_state
    use arg_parse

contains
    ! prepare a separate state with given dimensions
    !
    ! Inputs:
    !   N (integer): Number of subsystems
    !   D (integer): Dimension of each subsystem
    !   debug_level (integer): Debug level
    !
    ! Returns:
    !   state (complex*16 array): State vector
    !
    function prepare_separable_state(N, D, debug_level) result(state)
        implicit none

        integer N, D, dim
        real*8 norm
        integer ii, jj
        integer debug_level

        complex*16, dimension(:), allocatable :: state
        complex*16, dimension(:, :), allocatable :: separable_state
        integer, dimension(:), allocatable :: mat_indices

        dim = D ** N

        allocate(state(dim))
        allocate(separable_state(N, D))
        allocate(mat_indices(N))

        ! generate separable state containing N * d parameters,
        ! making sure to correctly normalize each wavefunction
        do ii = 1, N
            norm = 0
            do jj = 1, D
                separable_state(ii, jj) = cmplx(rand(0) * 2 - 1, rand(0) * 2 - 1)
                norm = norm + separable_state(ii, jj) * conjg(separable_state(ii, jj))
            end do
            separable_state(ii, :) = separable_state(ii, :) / sqrt(norm)
        end do

        if (debug_level .ge. DEBUG_LEVEL_1) then
            print *, "Wavefunction factorization = "
            call print_complex_matrix(separable_state)
        end if

        ! populate tensor
        state = 1

        do ii = 1, dim
            mat_indices = tensor2mat(ii, N, D)

            if (debug_level .ge. DEBUG_LEVEL_2) then
                print *, "Matrix indices for tensor index = ", ii, " are ", mat_indices
            end if

            do jj = 1, N
                state(ii) = state(ii) * separable_state(jj, mat_indices(jj))
            end do
        end do

        deallocate(separable_state, mat_indices)

    end function

    ! prepare a generic state with given dimensions
    !
    ! Inputs:
    !   N (integer): Number of subsystems
    !   D (integer): Dimension of each subsystem
    !
    ! Returns:
    !   state (complex*16 array): State vector
    !
    function prepare_generic_state(N, D) result(state)
        implicit none

        integer N, D, dim
        real*8 norm
        integer ii

        complex*16, dimension(:), allocatable :: state

        dim = D ** N

        allocate(state(dim))

        norm = 0

        do ii = 1, dim
            state(ii) = cmplx(rand(0) * 2 - 1, rand(0) * 2 - 1)
            norm = norm + state(ii) * conjg(state(ii))
        end do

        state = state / sqrt(norm)

    end function

    ! prepare a bell state
    !
    ! Returns:
    !   state (complex*16 array): State vector
    !
    function prepare_bell_state() result(state)
        implicit none

        real*8, parameter :: c = 1 / sqrt(2d0)
        complex*16 state(4)

        ! prepare |psi> = 1/sqrt(2) * (|01> - |10>)
        state(1) = (0d0, 0d0)
        state(2) = (c, 0d0)
        state(3) = -1.0 * (c, 0d0)
        state(4) = (0d0, 0d0)

    end function

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
    !   debug_level (integer): Debug level
    !
    ! Returns:
    !   rho_reduced_L (complex*16 matrix): Left reduced density matrix with dimensions (D ** (N - M), D ** (N - M))
    !
    function compute_left_reduced_density_matrix(rho, D, N, M, debug_level) result(rho_reduced_L)
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
                    if (debug_level .ge. DEBUG_LEVEL_2) then
                        print "('Added rho index ', (i4), (i4), ' to rho reduced L index ', (i4), (i4))", &
                            idx1, idx2, ii, jj
                    end if
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
    !   debug_level (integer): Debug level
    !
    ! Returns:
    !   rho_reduced_R (complex*16 matrix): Right reduced density matrix with dimensions (D ** M, D ** M)
    !
    function compute_right_reduced_density_matrix(rho, D, N, M, debug_level) result(rho_reduced_R)
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
                    if (debug_level .ge. DEBUG_LEVEL_2) then
                        print "('Added rho index ', (i4), (i4), ' to rho reduced R index ', (i4), (i4))", &
                            idx1, idx2, ii, jj
                    end if
                end do
            end do
        end do

    end function

    ! convert a single-dimensional tensor index into matrix form
    !
    ! Example: For a N=2 qubit system, you have
    !   tt = 1 => mat_idx = (1, 1)
    !   tt = 2 => mat_idx = (2, 1)
    !   tt = 3 => mat_idx = (1, 2)
    !   tt = 4 => mat_idx = (2, 2)
    !
    ! Inputs:
    !   tt (integer): Single-dimensional tensor index
    !   N (integer): Number of subsystems
    !   D (integer): Dimension of each subsystem
    !
    ! Returns:
    !   mat_idx (integer array): Array of indices into each subsystem
    !
    function tensor2mat(tt, N, D) result(mat_idx)
        implicit none

        integer tt, N, D
        integer ii
        integer mat_idx(N)

        do ii = 1, N
            mat_idx(ii) = modulo((tt - 1) / D ** (ii - 1), D) + 1
        end do

    end function

    ! get the norm of the given state vector
    !
    ! Inputs:
    !   state (complex*16 array): State to get norm of
    !
    ! Returns:
    !   norm (real*8): Norm
    !
    function get_norm(state) result(norm)
        implicit none

        integer ii
        real*8 norm
        complex*16, dimension(:) :: state

        norm = 0
        do ii = 1, size(state)
            norm = norm + state(ii) * conjg(state(ii))
        end do
    end function

    ! get the trace of the given density matrix
    !
    ! Inputs:
    !   rho (complex*16 array): Density matrix
    !
    ! Returns:
    !   tr (complex*16): Trace
    !
    function get_trace(rho) result(tr)
        implicit none

        integer ii
        complex*16 tr
        complex*16, dimension(:, :) :: rho

        tr = (0d0, 0d0)

        do ii = 1, size(rho, 1)
            tr = tr + rho(ii, ii)
        end do
    end function

    ! print complex vector in a nice format
    !
    ! Inputs:
    !   M (complex*16 array): Vector to print
    !
    subroutine print_complex_vector(M)
        implicit none

        integer ii
        complex*16, dimension(:) :: M

        do ii = 1, size(M)
            print '(*(sp, f7.4, 1x, f7.4, "i", 3x))', M(ii)
        end do
    end subroutine

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


program density_matrix
    use arg_parse
    use entropy
    use many_body_quantum_state
    implicit none

    real*8 S_L, S_R
    complex*16, dimension(:), allocatable :: state
    complex*16, dimension(:, :), allocatable :: rho
    complex*16, dimension(:, :), allocatable :: rho_reduced_L, rho_reduced_R

    ! variables to clock algorithm
    real*8 start, finish

    ! read x boundaries and number of discretization points
    call parse_cmd_args()
    write (arg_char, "(i8)") N
    print *, "N = ", adjustl(arg_char)
    write (arg_char, "(i8)") D
    print *, "D = ", adjustl(arg_char)
    print *, "type = ", system_type
    write (arg_char, "(i8)") M
    print *, "M = ", adjustl(arg_char)
    write (arg_char, "(i2)") debug_level
    print *, "debug = ", adjustl(arg_char)
    print *, "output_filename = ", output_filename

    ! prepare state
    if (system_type .eq. "separable") then
        state = prepare_separable_state(N, D, debug_level)
    else if (system_type .eq. "generic") then
        state = prepare_generic_state(N, D)
    else if (system_type .eq. "bell") then
        state = prepare_bell_state()
    end if

    ! compute density matrix
    rho = compute_density_matrix(state)

    ! compute left and right reduced density matrices
    rho_reduced_L = compute_left_reduced_density_matrix(rho, D, N, M, debug_level)
    rho_reduced_R = compute_right_reduced_density_matrix(rho, D, N, M, debug_level)

    ! compute entropy
    S_L = compute_entropy(rho_reduced_L, debug_level)
    S_R = compute_entropy(rho_reduced_R, debug_level)

    if (debug_level .ge. DEBUG_LEVEL_1) then
        ! print state
        print *, "State = "
        call print_complex_vector(state)
        print "('Norm = ', f6.4)", get_norm(state)

        ! print density matrices
        print *, "Density matrix = "
        call print_complex_matrix(rho)
    end if

    print "('Trace of density matrix = ', f6.4, 1x, sp, f7.4, 'i')", get_trace(rho)

    if (debug_level .ge. DEBUG_LEVEL_1) then
        ! print reduced density matrices
        print *, "Left reduced density matrix = "
        call print_complex_matrix(rho_reduced_L)

        print *, "Right reduced density matrix = "
        call print_complex_matrix(rho_reduced_R)
    end if

    print "('Trace of left reduced density matrix = ', f6.4, 1x, sp, f7.4, 'i')", get_trace(rho_reduced_L)
    print "('Trace of right reduced density matrix = ', f6.4, 1x, sp, f7.4, 'i')", get_trace(rho_reduced_R)

    ! print entropy
    print "('Entropy of left partition = ', f7.5)", S_L
    print "('Entropy of right partition = ', f7.5)", S_R

    ! write solution to file
    open(1, file=output_filename)
    write(1, *) "Entropy of left partition = ", S_L
    write(1, *) "Entropy of right partition = ", S_R
    write(1, *) "Dim density matrix = ", size(rho, 1)
    write(1, *) "Density matrix = "
    do ii = 1, size(rho, 1)
        write(1, *) rho(ii, :)
    end do
    write(1, *) "Dim left density matrix = ", size(rho_reduced_L, 1)
    write(1, *) "Left reduced density matrix = "
    do ii = 1, size(rho_reduced_L, 1)
        write(1, *) rho_reduced_L(ii, :)
    end do
    write(1, *) "Dim right density matrix = ", size(rho_reduced_R, 1)
    write(1, *) "Right reduced density matrix = "
    do ii = 1, size(rho_reduced_R, 1)
        write(1, *) rho_reduced_R(ii, :)
    end do
    close(1)

    deallocate(state, rho, rho_reduced_L, rho_reduced_R)

end program
