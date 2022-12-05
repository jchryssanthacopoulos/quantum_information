! Solve Split Operator
! ====================
!
! This program solves the time-dependent quantum harmonic oscillator using the split operator method
!
! Command-line arguments:
!   xmin (float): Minimum x value in the domain
!   xmax (float): Maximum x value in the domain
!   tmax (float): Maximum t value to simulate
!   num_x_pts (int): Number of points to discretize x coordinate
!   num_t_pts (int): Number of points to discretize t coordinate
!   output_filename (str): Name of file to save solution
!   debug (logical): Whether to print debug information
!
! Returns:
!   Saves the wavefunction and discretization parameters to the specified file
!


module split_operator

    real*8, parameter :: pi = 4.d0 * datan(1.d0)

contains
    ! evolve the state forward using the split operator method
    !
    ! Inputs:
    !   init_state (complex*16 array): Initial state
    !   final_state (complex*16 array | out): Final state
    !   x_grid (real*8 array): Grid of x values
    !   time (real*8): Time of initial state
    !   tmax (real*8): Maximum time
    !   Nx (integer): Number of x discretizaton points
    !   Nt (integer): Number of t discretization points
    !
    subroutine evolve_state(init_state, final_state, x_grid, time, tmax, Nx, Nt, debug)
        implicit none

        integer ii
        real*8 time
        logical debug

        ! discretization parameters
        integer Nx, Nt
        real*8 xmin, xmax, tmax
        real*8 dx, dp, dt
        real*8 x_grid(Nx)

        ! state variables
        complex*16, dimension(:), intent(in) :: init_state
        complex*16, dimension(:), intent(inout) :: final_state

        ! FFT variables
        integer*8 plan
        complex*16, dimension(:), allocatable :: state_transform

        ! to store kinetic and potential values
        real*8, dimension(:), allocatable :: T
        real*8, dimension(:), allocatable :: V

        xmin = x_grid(1)
        xmax = x_grid(Nx)

        allocate(state_transform(Nx))
        allocate(T(Nx))
        allocate(V(Nx))

        ! get spacing in space, momentum, and time
        dx = (xmax - xmin) / (Nx - 1)
        dp = 2 * pi / (xmax - xmin)
        dt = tmax / (Nt - 1)

        ! multiply by potential part of Hamiltonian
        do ii = 1, Nx
            V(ii) = potential(x_grid(ii), time, tmax)
            final_state(ii) = cexp(cmplx(0.0, -0.5 * V(ii) * dt)) * init_state(ii)
        end do

        ! normalize
        call normalize(final_state, dx)

        ! call FFT to go from coordinate space to momentum space
        call dfftw_plan_dft_1d(plan, Nx, final_state, state_transform, -1, 64)
        call dfftw_execute_dft(plan, final_state, state_transform)
        call dfftw_destroy_plan(plan)

        ! make sure to normalize in momentum space
        final_state = state_transform
        call normalize(final_state, dp)

        ! multiply by kinetic part of Hamiltonian
        do ii = 1, Nx
            T(ii) = kinetic(ii, Nx, xmin, xmax)
            final_state(ii) = cexp(cmplx(0.0, -1.0 * T(ii) * dt)) * final_state(ii)
        end do

        call normalize(final_state, dp)

        ! call inverse FFT to go from momentum space to coordinate space
        call dfftw_plan_dft_1d(plan, Nx, final_state, state_transform, 1, 64)
        call dfftw_execute_dft(plan, final_state, state_transform)
        call dfftw_destroy_plan(plan)

        final_state = state_transform
        call normalize(final_state, dx)

        ! multiply by potential part of Hamiltonian
        do ii = 1, Nx
            final_state(ii) = cexp(cmplx(0.0, -0.5 * V(ii) * dt)) * final_state(ii)
        end do

        call normalize(final_state, dx)

        if (debug) then
            print *, "time = ", time
            print *, "T = ", T
            print *, "V = ", V
        end if

        deallocate(state_transform, T, V)

    end subroutine

    ! get potential for given position and time
    !
    ! Inputs:
    !   q (real*8): Position
    !   t (real*8): Time
    !   tmax (real*8): Maximum time
    !
    ! Returns:
    !   V (real*8): Potential at given position and time
    !
    function potential(q, t, tmax) result(V)
        implicit none

        real*8 q, t, tmax
        real*8 V

        V = 0.5 * (q - t / tmax) ** 2
    end function

    ! get kinetic term at given bin value
    !
    ! Inputs:
    !   bin (integer): Bin
    !   Nx (integer): Number of total bins
    !   xmin (real*8): Minimum x value
    !   xmax (real*8): Maximum x value
    !
    ! Returns:
    !   T (real*8): Kinetic term
    !
    function kinetic(bin, Nx, xmin, xmax) result(T)
        implicit none

        integer bin, Nx
        real*8 xmin, xmax, p, T

        if (bin <= int(Nx / 2)) then
            p = 2 * pi * (bin - 1) / (xmax - xmin)
        else if (bin > int(Nx / 2)) then
            p = 2 * (pi * (bin - Nx - 1)) / (xmax - xmin)
        end if

        T = 0.5 * p ** 2

    end function

    ! normalize function using integral with given step
    !
    ! Inputs:
    !   psi: Function to normalize
    !   step: Step size in the domain to integrate over
    !
    subroutine normalize(psi, step)
        implicit none

        complex*16, dimension(:), intent(inout) :: psi
        real*8 step, norm

        norm = sum(psi * conjg(psi) * step)
        psi = psi / sqrt(norm)

    end subroutine

end module


program solve_time_dep_ho
    use arg_parse
    use split_operator
    implicit none

    ! discretization parameters
    real*8 dx, dt
    real*8, dimension(:), allocatable :: x_grid, t_grid

    ! to save wavefunction as a function of time
    complex*16, dimension(:,:), allocatable :: psi

    ! read x boundaries and number of discretization points
    call parse_cmd_args()
    write (arg_char, "(f7.3)") xmin
    print *, "xmin = ", adjustl(arg_char)
    write (arg_char, "(f7.3)") xmax
    print *, "xmax = ", adjustl(arg_char)
    write (arg_char, "(f7.3)") tmax
    print *, "tmax = ", adjustl(arg_char)
    write (arg_char, "(i8)") num_x_pts
    print *, "num_x_pts = ", adjustl(arg_char)
    write (arg_char, "(i8)") num_t_pts
    print *, "num_t_pts = ", adjustl(arg_char)
    write (arg_char, "(l1)") debug
    print *, "debug = ", adjustl(arg_char)
    print *, "output_filename = ", output_filename

    allocate(x_grid(num_x_pts))
    allocate(t_grid(num_t_pts))
    allocate(psi(num_x_pts, num_t_pts))

    ! get spacing
    dx = (xmax - xmin) / (num_x_pts - 1)
    dt = tmax / (num_t_pts - 1)

    psi = 0

    ! initialize state with ground state of harmonic oscillator
    do ii = 1, num_x_pts
        x_grid(ii) = xmin + (ii - 1) * dx
        psi(ii, 1) = pi ** (-0.25d0) * exp(-x_grid(ii) ** 2d0 / 2d0)
    end do

    ! normalize
    call normalize(psi(:, 1), dx)

    ! generate t lattice
    do ii = 1, num_t_pts
        t_grid(ii) = (ii - 1) * dt
    end do

    if (debug) then
        print *, "x_grid = ", x_grid
        print *, "t_grid = ", t_grid
    end if

    ! propagate state
    do ii = 2, num_t_pts
        call evolve_state(psi(:, ii - 1), psi(:, ii), x_grid, t_grid(ii), tmax, num_x_pts, num_t_pts, debug)
    end do

    ! write solution to file
    open(1, file=output_filename)
    write(1, *) "x grid =", x_grid
    write(1, *) "t grid =", t_grid
    do ii = 1, num_t_pts
        write(1, *) psi(:, ii)
    end do
    close(1)

    deallocate(x_grid, t_grid, psi)

end program
