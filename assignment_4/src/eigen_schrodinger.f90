! Eigen Schrodinger
! =================
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
    allocate(x_grid(npoints + 1))
    allocate(H_diag(npoints + 1))
    allocate(H_off_diag(npoints))
    allocate(eigenvectors(npoints + 1, npoints + 1))
    allocate(work(2 * npoints))

    ! get spacing
    dx = (xmax - xmin) / npoints

    H_off_diag = -1 / (2 * dx ** 2)

    do ii = 1, npoints + 1
        x_grid(ii) = xmin + (ii - 1) * dx
        H_diag(ii) = 1 / (dx ** 2) + 0.5 * x_grid(ii) ** 2
    end do

    ! compute eigenvalues and eigenvectors
    call dsteqr("I", npoints + 1, H_diag, H_off_diag, eigenvectors, npoints + 1, work, info)

    if (info .eq. 0) then
        print *, "Success!"
    else
        print *, "Failed to obtain eigenvalues and eigenvectors"
    end if

    ! write solution to file
    open(1, file=output_filename)
    write(1, *) "N =", npoints + 1
    write(1, *) "x grid =", x_grid
    do ii = 1, npoints + 1
        write(1, *) "Eigenvalue =", H_diag(ii)
        write(1, *) "Eigenvector =" 
        write(1, *) eigenvectors(ii, :)
    end do
    close(1)

    ! free memory
    deallocate(x_grid, H_diag, H_off_diag, eigenvectors, work)

end program
