! Eigen Schrodinger
! =================
!
! This program solves the time-independent Schrodinger equation for a harmonic oscillator through discretization
!
! Command-line arguments:
!   xmin (float): Minimum x value in the domain
!   xmax (float): Maximum x value in the domain
!   npoints: Number of points to discretize domain
!   output_filename: Name of file to save eigenvalues and eigenvectors
!
! Returns:
!   Saves number of points, x grid, eigenvalues, and eigenvectors to the specified file
!


program eigen_schrodinger
    use arg_parse
    implicit none

    ! discretization parameters
    real*8 dx
    real*8, dimension(:), allocatable :: x_grid

    ! variables to get eigenvalues and eigenvectors
    integer info
    real*8, dimension(:), allocatable :: H_diag, H_off_diag
    real*8, dimension(:,:), allocatable :: eigenvectors
    real*8, dimension(:), allocatable :: work

    ! read x boundaries and number of discretization points
    call parse_cmd_args()
    write (arg_char, "(f7.3)") xmin
    print *, "xmin = ", adjustl(arg_char)
    write (arg_char, "(f7.3)") xmax
    print *, "xmax = ", adjustl(arg_char)
    write (arg_char, "(i8)") npoints
    print *, "npoints = ", adjustl(arg_char)
    print *, "output_filename = ", output_filename

    ! allocate memory
    allocate(x_grid(npoints))
    allocate(H_diag(npoints))
    allocate(H_off_diag(npoints - 1))
    allocate(eigenvectors(npoints, npoints))
    allocate(work(2 * (npoints - 1)))

    ! get spacing
    dx = (xmax - xmin) / (npoints - 1)

    ! discretize Hamiltonian
    H_off_diag = -1 / (2 * dx ** 2)

    do ii = 1, npoints
        x_grid(ii) = xmin + (ii - 1) * dx
        H_diag(ii) = 1 / (dx ** 2) + 0.5 * x_grid(ii) ** 2
    end do

    ! compute eigenvalues and eigenvectors
    call dsteqr("I", npoints, H_diag, H_off_diag, eigenvectors, npoints, work, info)

    if (info .eq. 0) then
        print *, "Success!"
    else
        print *, "Failed to obtain eigenvalues and eigenvectors"
    end if

    ! write solution to file
    open(1, file=output_filename)
    write(1, *) "N =", npoints
    write(1, *) "x grid =", x_grid
    do ii = 1, npoints
        write(1, *) "Eigenvalue =", H_diag(ii)
        write(1, *) "Eigenvector =" 
        write(1, *) eigenvectors(ii, :)
    end do
    close(1)

    ! free memory
    deallocate(x_grid, H_diag, H_off_diag, eigenvectors, work)

end program
