module cmatrix_type
    implicit none

    type cmatrix
        integer, dimension(2) :: N
        complex*16, dimension(:, :), allocatable :: elems
        complex*16 trace
    end type

    ! this is not an operator since we want it to modify the trace of the type itself
    interface Trace
        module procedure Trace
    end interface

    interface operator (.Adj.)
        module procedure Adjoint
    end interface

contains
    ! initialize the matrix with random data
    function Init(nrows, ncols) result(M)
        integer*4 nrows, ncols
        type(cmatrix) M

        ! matrices to store real and complex parts of initialized matrix
        real*8, dimension(:, :), allocatable :: a, b

        ! allocate matrices
        allocate(M%elems(nrows, ncols))
        allocate(a(nrows, ncols))
        allocate(b(nrows, ncols))

        ! initialize real and complex parts with random data
        call random_number(a)
        call random_number(b)

        M%N = (/nrows, ncols/)
        M%elems = cmplx(a, b, kind(1d0))

        ! calculate trace
        call Trace(M)
    end function

    ! calculate the trace and store in the variable "trace"
    subroutine Trace(M)
        implicit none

        type(cmatrix), intent(inout) :: M
        integer*4 ii

        ! initialize trace to zero
        M%trace = (0d0, 0d0)

        if (M%N(1) /= M%N(2)) then
            print *, "Cannot compute trace for non-square matrix"
        end if

        do ii = 1, M%N(1)
            M%trace = M%trace + M%elems(ii, ii)
        end do
    end subroutine

    ! calculate the matrix adjoint and return in a new variable
    function Adjoint(M) result(Mdagger)
        type(cmatrix), intent(in) :: M
        type(cmatrix) Mdagger

        ! set transposed dimensions
        Mdagger%N = (/M%N(2), M%N(1)/)

        ! take complex conjugate transpose
        Mdagger%elems = conjg(transpose(M%elems))

        ! the trace is given by the complex conjugate
        Mdagger%trace = conjg(M%trace)
    end function

    ! deallocate memory associated with matrix elements
    subroutine Del(M)
        implicit none

        type(cmatrix) M

        if (Allocated(M%elems)) then
            deallocate(M%elems)
            M%N = 0  ! set dimensions to zero
        end if
    end subroutine

    ! print the matrix in a nice format
    subroutine print_matrix(M)
        type(cmatrix) M
        integer*4 ii

        do ii = 1, M%N(1)
            print '(*(sp, f8.4, 1x, f7.4, "i", 3x))', M%elems(ii, :)
        end do
    end subroutine

    ! print the trace of the matrix only if it's been calculated
    subroutine print_trace(M)
        type(cmatrix) M

        if (M%N(1) /= M%N(2)) then
            print *, "Cannot display trace for non-square matrix"
            return
        end if

        print '(1x, "The trace of M is ", sp, f7.4, 1x, f7.4, "i")', M%trace
    end subroutine

    ! write matrix to file
    subroutine write_matrix(M, filename)
        implicit none

        type(cmatrix) M
        integer*4 ii
        character(*) filename

        open(1, file=filename)

        ! write dimensions
        write(1, '(A, i5.1, A, i5.1, /)') 'Dimensions: ', M%N(1), ' x ', M%N(2)

        ! write matrix elements
        write(1, *) "Matrix elements:"
        do ii = 1, M%N(1)
            write(1, '(*(sp, f8.4, 1x, f7.4, "i", 2x))')  M%elems(ii, :)
        end do

        ! write trace
        write(1, *) "Trace:"
        if (M%N(1) /= M%N(2)) then
            write(1, *) "Could not calculate trace for non-square matrix"
        else
            write(1, '(sp, f8.4, 1x, f7.4, "i")') M%trace
        end if

        close(1)
    end subroutine

end module
