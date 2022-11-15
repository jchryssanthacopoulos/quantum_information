! Exercise 2
! ==========
!


module histogram

contains
    subroutine make_hist(values, nbins, min_val, max_val, filename)
        implicit none

        integer nbins
        real*8 min_val, max_val
        real*8, dimension(:), intent(in) :: values
        character(len=*), intent(in) :: filename

        integer ii, jj
        integer bin_count
        real*8 dx, norm_factor, left_bin_edge

        ! compute spacing and normalization factor
        dx = (max_val - min_val) / nbins
        norm_factor = size(values) * dx

        open(1, file=filename)
        write(1, *) "left_edges centers count norm_count"

        ! compute histogram
        left_bin_edge = min_val

        do ii = 1, nbins
            bin_count = 0

            do jj = 1, size(values, 1)
                if (values(jj) >= left_bin_edge .and. values(jj) <= left_bin_edge + dx) then
                    bin_count = bin_count + 1
                end if
            end do

            ! write entry to file
            write(1, *) left_bin_edge, left_bin_edge + dx, bin_count, bin_count / norm_factor

            left_bin_edge = left_bin_edge + dx
        end do

        close(1)
    end subroutine

end module


module arg_parse
    implicit none

    integer ii

    integer ndim, nsamples, nbins
    real*8 min_val, max_val
    logical debug_mode

    ! default paramters
    integer, parameter :: ndim_default = 10
    integer, parameter :: nsamples_default = 1
    integer, parameter :: nbins_default = 100
    real*8, parameter :: min_val_default = 0.0
    real*8, parameter :: max_val_default = 5.0
    logical, parameter :: debug_mode_default = .false.

contains
    ! parse command-line arguments
    subroutine parse_cmd_args
        implicit none

        integer num_args
        character(len=32) arg

        ! set defaults
        ndim = ndim_default
        nsamples = nsamples_default
        nbins = nbins_default
        min_val = min_val_default
        max_val = max_val_default
        debug_mode = debug_mode_default

        num_args = command_argument_count()

        ii = 0

        do while(ii <= num_args)
            call get_command_argument(ii, arg)

            select case (arg)
                case ('--ndim')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) ndim
                        ii = ii + 1
                    end if

                case ('--nsamples')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) nsamples
                        ii = ii + 1
                    end if

                case ('--nbins')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) nbins
                        ii = ii + 1
                    end if

                case ('--min_val')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) min_val
                        ii = ii + 1
                    end if

                case ('--max_val')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) max_val
                        ii = ii + 1
                    end if

                case ('-d', '--debug')
                    debug_mode = .true.
            end select

            ii = ii + 1
        end do

    end subroutine

end module


program exercise_2
    use histogram
    use arg_parse
    implicit none

    ! to store hermitian matrix
    complex*16, dimension(:, :), allocatable :: H

    ! variables for computing eigenvalues and derived quantities
    integer lwork, info
    real*8 ave_delta_eigvals
    real*8, dimension(:), allocatable :: eigvals, norm_eigval_spacings
    complex*16, dimension(:), allocatable :: work
    real*8, dimension(:), allocatable :: rwork

    ! variables for storing eigenvalue spacings across samples
    integer ss
    integer min_samp_range, max_samp_range
    real*8, dimension(:), allocatable :: all_norm_eigval_spacings

    ! read matrix dimension and histogram parameters
    call parse_cmd_args()
    print *, "ndim =", ndim
    print *, "nsamples =", nsamples
    print *, "nbins =", nbins
    print *, "min_val =", min_val
    print *, "max_val =", max_val

    ! used for computing eigenvalues
    lwork = max(1, 2 * ndim - 1)

    ! allocate memory
    allocate(H(ndim, ndim))
    allocate(eigvals(ndim))
    allocate(norm_eigval_spacings(ndim - 1))
    allocate(all_norm_eigval_spacings(nsamples * (ndim - 1)))
    allocate(work(max(1, lwork)))
    allocate(rwork(max(1, 3 * ndim - 2)))

    do ss = 1, nsamples
        ! generatre random hermitian matrix
        H = rand_hermitian_matrix(ndim)

        ! display matrix
        if (debug_mode) then
            print *, "Sampled Hermitian matrix ="
            call print_matrix(H, ndim)
        end if

        ! compute eigenvalues in ascending order
        call zheev('N', 'U', ndim, H, ndim, eigvals, work, lwork, rwork, info)

        ! compute eigenvalue spacings and average
        ave_delta_eigvals = (eigvals(ndim) - eigvals(1)) / (ndim - 1)
        do ii = 1, ndim - 1
            norm_eigval_spacings(ii) = (eigvals(ii + 1) - eigvals(ii)) / ave_delta_eigvals
        end do

        ! print the eigenvalues
        if (info .ne. 0) then
            print *, "Encountered error code when computing eigenvalues: ", info
            return
        end if

        if (debug_mode) then
            print *, "Eigenvalues =", eigvals
            print *, "Normalized eigenvalue spacings =", norm_eigval_spacings
            print *, "Average eigenvalue spacing =", ave_delta_eigvals
        end if

        ! add to total samples
        min_samp_range = 1 + (ss - 1) * (ndim - 1)
        max_samp_range = ss * (ndim - 1)
        all_norm_eigval_spacings(min_samp_range:max_samp_range) = norm_eigval_spacings
    end do

    ! compute histogram
    call make_hist(all_norm_eigval_spacings, nbins, min_val, max_val, "histogram.txt")

    deallocate(H, eigvals, norm_eigval_spacings, all_norm_eigval_spacings, work, rwork)

contains
    function rand_hermitian_matrix(n) result(H)
        integer n
        integer ii, jj
        complex*16, dimension(n, n) :: H

        do ii = 1, n
            do jj = 1, ii
                if (ii /= jj) then
                    ! sample random complex numbers off the diagonal
                    h(ii, jj) = cmplx(RAND(0)*2 - 1, RAND(0)*2 - 1)
                    h(jj, ii) = conjg(h(ii, jj))
                else
                    ! sample real numbers on the diagonal
                    h(ii, ii) = RAND(0)*2 - 1
                end if
            end do
        end do
    end function

    subroutine print_matrix(M, n)
        implicit none

        integer n
        integer ii
        complex*16 M(n, n)

        do ii = 1, n
            print '(*(sp, f7.4, 1x, f7.4, "i", 3x))', M(ii, :)
        end do
    end subroutine

end program
