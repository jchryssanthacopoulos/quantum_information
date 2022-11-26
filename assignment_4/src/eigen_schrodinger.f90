! Eigen Schrodinger
! =================
!


program eigen_schrodinger
    implicit none

    integer ii

    ! discretization parameters
    integer N
    real*8 x_min, x_max, dx
    real*8, dimension(:), allocatable :: x_grid

    ! Hamiltonian parameters
    real*8 omega

    ! variables to get eigenvalues and eigenvectors
    integer info
    real*8, dimension(:), allocatable :: H_diag, H_off_diag
    real*8, dimension(:,:), allocatable :: eigenvectors
    real*8, dimension(:), allocatable :: work

    ! set Hamiltonian parameters
    omega = 1.0

    ! number of points
    N = 1000
    x_min = -5.0
    x_max = 5.0

    ! allocate memory
    allocate(x_grid(N + 1))
    allocate(H_diag(N + 1))
    allocate(H_off_diag(N))
    allocate(eigenvectors(N + 1, N + 1))
    allocate(work(2 * N))

    ! get spacing
    dx = (x_max - x_min) / N

    H_off_diag = -1 / (2 * dx ** 2)

    do ii = 1, N + 1
        x_grid(ii) = x_min + (ii - 1) * dx
        H_diag(ii) = 1 / (dx ** 2) + 0.5 * (omega * x_grid(ii)) ** 2
    end do

    ! compute eigenvalues and eigenvectors
    call dsteqr("I", N + 1, H_diag, H_off_diag, eigenvectors, N + 1, work, info)

    print *, "info =", info

    ! normalize the eigenvectors
    ! do ii = 1, N + 1
    !     eigenvectors(:, ii) = eigenvectors(:, ii) / sqrt(sum(eigenvectors(:, ii) ** 2))
    ! end do

    ! write grid and eigenvectors to file
    open(1, file="eigenvectors.txt")
    write(1, *) x_grid
    do ii = 1, N + 1
        write(1, *) eigenvectors(ii, :)
    end do
    close(1)

    ! write eigenvalues to file
    open(1, file="eigenvalues.txt")
    write(1, *) H_diag
    close(1)

    ! free memory
    deallocate(x_grid, H_diag, H_off_diag, eigenvectors, work)

end program
