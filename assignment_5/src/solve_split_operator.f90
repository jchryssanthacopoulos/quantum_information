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
!


program solve_split_operator
    use arg_parse
    implicit none

    ! constants
    real*8 pi

    ! discretization parameters
    real*8 dx, dt
    real*8, dimension(:), allocatable :: x_grid, t_grid

    ! FFT variables
    integer*8 plan
    complex*8, dimension(:), allocatable :: psi_transform_in
    complex*8, dimension(:), allocatable :: psi_transform_out

    ! to save wavefunction as a function of time
    complex*8, dimension(:,:), allocatable :: psixt

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
    print *, "output_filename = ", output_filename

    allocate(x_grid(num_x_pts))
    allocate(t_grid(num_t_pts))
    allocate(psixt(num_x_pts, num_t_pts))
    allocate(psi_transform_in(num_x_pts))
    allocate(psi_transform_out(num_x_pts))

    pi = 4.d0 * datan(1.d0)

    ! get spacing
    dx = (xmax - xmin) / (num_x_pts - 1)
    dt = tmax / (num_t_pts - 1)

    psixt = 0
    psi_transform_in = 0
    psi_transform_out = 0

    ! initialize state with ground state of harmonic oscillator
    do ii = 1, num_x_pts
        x_grid(ii) = xmin + (ii - 1) * dx
        psixt(ii, 1) = pi ** (-0.25d0) * exp(-x_grid(ii) ** 2d0 / 2d0)
    end do

    ! call FFT on ground state
    psi_transform_in = psixt(:, 1)
    call dfftw_plan_dft_1d(plan, num_x_pts, psi_transform_in, psi_transform_out, -1, 64)
    call dfftw_execute_dft(plan, psi_transform_in, psi_transform_out)
    call dfftw_destroy_plan(plan)
    psixt(:, 2) = psi_transform_out

    ! call inverse FFT
    psi_transform_in = psi_transform_out
    psi_transform_out = 0
    call dfftw_plan_dft_1d(plan, num_x_pts, psi_transform_in, psi_transform_out, 1, 64)
    call dfftw_execute_dft(plan, psi_transform_in, psi_transform_out)
    call dfftw_destroy_plan(plan)
    psixt(:, 3) = psi_transform_out

    do ii = 1, num_t_pts
        t_grid(ii) = (ii - 1) * dt
    end do

    ! write solution to file
    open(1, file=output_filename)
    write(1, *) "x grid =", x_grid
    write(1, *) "t grid =", t_grid
    do ii = 1, num_t_pts
        write(1, *) psixt(:, ii)
    end do
    close(1)

    deallocate(x_grid, t_grid, psixt, psi_transform_in, psi_transform_out)

end program
