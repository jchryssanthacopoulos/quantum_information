! Exercise 2
! ==========
!
! This program computes the normalized eigenvalue spacings for a given matrix type and saves them to a file
!
! Command-line arguments:
!   mat_type: Type of matrix to diagonalize (options are `hermitian` and `diag`) (default = "hermitian")
!   output_filename: Name of file to save the output histogram (default = "histogram.csv")
!   ndim: Number of rows and columns of the matrix (default = 10)
!   nsamples: Number of random matrix to sample to produce histogram (default = 1)
!   nbins: Number of bins to use to produce histogram (default = 100)
!   min_val: Minimum spacings value of histogram (default = 0.0)
!   max_val: Maximum spacings value of histogram (default = 5.0)
!   debug: Flag indicating whether to run in debug mode, which prints the matrix, eigenvalues, spacings, and
!     their average for each sample (default = False)
!
! Returns:
!   File of normalized eigenvalue spacings
!   If debug flag passed, matrices, eigenvalues, and other information are also displayed
!


! this module parses the command-line arguments needed to compute the histogram
module arg_parse_eigenvalues
    implicit none

    integer ii
    
    character(len=50) mat_type
    character(len=50) output_filename
    integer ndim, nsamples, nbins
    real*8 min_val, max_val
    logical debug_mode

    ! default parameters
    character(len=50), parameter :: mat_type_default = "hermitian"
    character(len=50), parameter :: output_filename_default = "histogram.csv"
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
        mat_type = mat_type_default
        output_filename = output_filename_default
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
                case ('--mat_type')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, *) mat_type
                        ii = ii + 1
                    end if

                case ('--output_filename')
                    if (ii < num_args) then
                        call get_command_argument(ii + 1, arg)
                        read(arg, '(A)') output_filename
                        ii = ii + 1
                    end if

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
    use arg_parse_eigenvalues
    use histogram
    use mat_ops
    implicit none

    ! for displaying command-line arguments
    character(len=8) arg_char

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
    print *, "mat_type = ", mat_type
    print *, "output_filename = ", output_filename
    write (arg_char, "(i8)") ndim
    print *, "ndim = ", adjustl(arg_char)
    write (arg_char, "(i8)") nsamples
    print *, "nsamples = ", adjustl(arg_char)
    write (arg_char, "(i8)") nbins
    print *, "nbins = ", adjustl(arg_char)
    write (arg_char, "(f6.3)") min_val
    print *, "min_val = ", adjustl(arg_char)
    write (arg_char, "(f6.3)") max_val
    print *, "max_val = ", adjustl(arg_char)

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
        ! generatre random hermitian or real diagonal matrix
        if (mat_type .eq. "hermitian") then
            H = rand_hermitian_matrix(ndim)
        else if (mat_type .eq. "diag") then
            H = rand_real_diag_matrix(ndim)
        end if

        ! display matrix
        if (debug_mode) then
            print *, "Sampled matrix ="
            call print_complex_matrix(H, ndim)
        end if

        ! compute eigenvalues in ascending order
        call zheev("N", "U", ndim, H, ndim, eigvals, work, lwork, rwork, info)

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
    call make_hist(all_norm_eigval_spacings, nbins, min_val, max_val, output_filename)

    deallocate(H, eigvals, norm_eigval_spacings, all_norm_eigval_spacings, work, rwork)

end program
