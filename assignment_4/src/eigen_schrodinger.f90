! Eigen Schrodinger
! =================
!


! this module parses the command-line arguments needed to discretize the eigenvalue problem
module arg_parse
    implicit none

    integer ii

    real*8 xmin, xmax
    integer npoints
    character(len=50) eigenvalues_filename
    character(len=50) eigenvectors_filename

    ! default parameters
    real*8, parameter :: xmin_default = -5.0
    real*8, parameter :: xmax_default = 5.0
    integer, parameter :: npoints_default = 1000
    character(len=50), parameter :: eigenvalues_filename_default = "eigenvalues.txt"
    character(len=50), parameter :: eigenvectors_filename_default = "eigenvectors.txt"

    ! for displaying command-line arguments
    character(len=8) arg_char

contains
    ! parse command-line arguments
    subroutine parse_cmd_args
        implicit none

        integer num_args
        character(len=32) arg

        ! set defaults
        xmin = xmin_default
        xmax = xmax_default
        npoints = npoints_default
        eigenvalues_filename = eigenvalues_filename_default
        eigenvectors_filename = eigenvectors_filename_default

        num_args = command_argument_count()

        ii = 0

        do while(ii <= num_args)
            call get_command_argument(ii, arg)

            select case (arg)
                case ("--xmin")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) xmin
                        ii = ii + 1
                    end if

                case ("--xmax")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) xmax
                        ii = ii + 1
                    end if

                case ("--npoints")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) npoints
                        ii = ii + 1
                    end if

                case ("--eigenvalues_filename")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, "(A)") eigenvalues_filename
                        ii = ii + 1
                    end if

                case ("--eigenvectors_filename")
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, "(A)") eigenvectors_filename
                        ii = ii + 1
                    end if

            end select

            ii = ii + 1
        end do

    end subroutine

end module


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
    write (arg_char, "(f6.3)") xmin
    print *, "xmin = ", adjustl(arg_char)
    write (arg_char, "(f6.3)") xmax
    print *, "xmax = ", adjustl(arg_char)
    write (arg_char, "(i8)") npoints
    print *, "npoints = ", adjustl(arg_char)
    print *, "eigenvalues_filename = ", eigenvalues_filename
    print *, "eigenvectors_filename = ", eigenvectors_filename

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

    print *, "info =", info

    ! normalize the eigenvectors
    ! do ii = 1, N + 1
    !     eigenvectors(:, ii) = eigenvectors(:, ii) / sqrt(sum(eigenvectors(:, ii) ** 2))
    ! end do

    ! write grid and eigenvectors to file
    open(1, file=eigenvectors_filename)
    write(1, *) x_grid
    do ii = 1, npoints + 1
        write(1, *) eigenvectors(ii, :)
    end do
    close(1)

    ! write eigenvalues to file
    open(1, file=eigenvalues_filename)
    write(1, *) H_diag
    close(1)

    ! free memory
    deallocate(x_grid, H_diag, H_off_diag, eigenvectors, work)

end program
