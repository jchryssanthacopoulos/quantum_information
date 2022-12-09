! Generate density matrix
! =======================
!
! This program generates the density matrix for a many-body quantum state and its reduced matrix
!
! Command-line arguments:
!   N (integer): Number of subsystems in the quantum state
!   dim (integer): Number of dimensions of each subsystem
!   type (character): Type of system to compute (options are "separable" and "generic")
!   M (integer): Number of subsystems to trace over to compute left and right reduced density matrices
!   output_filename (character): Name of file to save density matrices
!   debug (logical): Whether to print debug information
!
! Returns:
!   Saves the density matrices to the specified file
!


module many_body_quantum_state

contains

    ! generate a separate wavefunction state with given dimensions
    !
    ! Inputs:
    !   N (integer): Number of subsystems
    !   D (integer): Dimension of each subsystem
    !
    ! Returns:
    !   state (complex*16 array): Wavefunction vector
    !
    function gen_separable_state(N, D, debug) result(state)
        implicit none

        integer N, D, dim
        real*8 norm
        integer ii, jj
        logical debug

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

        ! populate tensor
        state = 1

        do ii = 1, dim
            mat_indices = tensor2mat(ii, N, D)

            if (debug) then
                print *, "Matrix indices for tensor index = ", ii, " are ", mat_indices
            end if

            do jj = 1, N
                state(ii) = state(ii) * separable_state(jj, mat_indices(jj))
            end do
        end do

        deallocate(separable_state, mat_indices)

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
    use many_body_quantum_state
    implicit none

    complex*16, dimension(:), allocatable :: state
    complex(8), dimension(:, :), allocatable :: rho

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
    write (arg_char, "(l1)") debug
    print *, "debug = ", adjustl(arg_char)
    print *, "output_filename = ", output_filename

    ! prepare state
    if (system_type .eq. "separable") then
        state = gen_separable_state(N, D, debug)
    end if

    ! compute density matrix
    rho = compute_density_matrix(state)

    if (debug) then
        ! print state
        print *, "State = "
        call print_complex_vector(state)
        print "('Norm = ', f6.4)", get_norm(state)

        ! print density matrix
        print *, "Density matrix = "
        call print_complex_matrix(rho)
        print "('Trace = ', f6.4, 1x, sp, f7.4, 'i')", get_trace(rho)
    end if

    deallocate(state, rho)

end program
