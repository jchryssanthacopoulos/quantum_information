! this module contains utilities related to performing operations on matrices
module mat_ops

contains
    ! multiply two matrices in row-column order
    !
    ! Inputs:
    !   A: First matrix to multiply
    !   B: Second matrix to multiply
    !   C: Matrix to store result
    !   dims: Number of rows of A, columns of B, and inner dimension
    !
    subroutine matmul_row_col(A, B, C, dims)
        implicit none

        integer*4 ii, jj, kk
        integer*4 dims(3)
        real*8 A(dims(1), dims(3)), B(dims(3), dims(2))
        real*8, intent(out) :: C(dims(1), dims(3))

        do ii = 1, dims(1)
            do jj = 1, dims(2)
                do kk = 1, dims(3)
                    C(ii, jj) = C(ii, jj) + A(ii, kk) * B(kk, jj)
                end do
            end do
        end do
    end subroutine

    ! multiply two matrices in column-row order
    !
    ! Inputs:
    !   A: First matrix to multiply
    !   B: Second matrix to multiply
    !   C: Matrix to store result
    !   dims: Number of rows of A, columns of B, and inner dimension
    !
    subroutine matmul_col_row(A, B, C, dims)
        implicit none

        integer*4 ii, jj, kk
        integer*4 dims(3)
        real*8 A(dims(1), dims(3)), B(dims(3), dims(2))
        real*8, intent(out) :: C(dims(1), dims(3))

        do kk = 1, dims(3)
            do jj = 1, dims(2)
                do ii = 1, dims(1)
                    C(ii, jj) = C(ii, jj) + A(ii, kk) * B(kk, jj)
                end do
            end do
        end do
    end subroutine

    ! display maximum absolute error between two matrices
    !
    ! Inputs:
    !   A: First matrix
    !   B: Second matrix
    !   nrows: Number of rows of matrices
    !   ncols: Number of columns
    !
    function max_abs_error(A, B, nrows, ncols)
        integer*4 nrows, ncols, ii, jj, cnt
        real*8 error, max_abs_error
        real*8 A(nrows, ncols), B(nrows, ncols)

        cnt = 0

        ! iterate over all entries, computing the absolute difference
        do ii = 1, nrows
            do jj = 1, ncols
                error = abs(A(ii, jj) - B(ii, jj))
                if (cnt == 0) then
                    max_abs_error = error
                    cnt = cnt + 1
                else if (error > max_abs_error) then
                    max_abs_error = error
                end if
            end do
        end do
    end function

    ! generate random hermitian matrix
    !
    ! Inputs:
    !   n: Number of rows and columns of matrix
    !
    ! Returns:
    !   H: Generated matrix
    !
    function rand_hermitian_matrix(n) result(H)
        integer n
        integer ii, jj
        complex*16 H(n, n)

        do ii = 1, n
            do jj = 1, ii
                if (ii /= jj) then
                    ! sample random complex numbers off the diagonal
                    H(ii, jj) = cmplx(rand(0)*2 - 1, rand(0)*2 - 1)
                    H(jj, ii) = conjg(H(ii, jj))
                else
                    ! sample real numbers on the diagonal
                    H(ii, ii) = rand(0)*2 - 1
                end if
            end do
        end do
    end function

    ! generate random real diagonal matrix
    !
    ! Inputs:
    !   n: Number of rows and columns of matrix
    !
    ! Returns:
    !   H: Generated matrix
    !
    function rand_real_diag_matrix(n) result(H)
        integer n
        integer ii
        complex*16, dimension(n, n) :: H

        H = 0

        do ii = 1, n
            ! sample real numbers on the diagonal
            H(ii, ii) = rand(0)*2 - 1
        end do
    end function

    ! print real matrix in a nice format
    !
    ! Inputs:
    !   M: Matrix to print
    !   nrows: Number of rows
    !   ncols: Number of columns
    !
    subroutine print_real_matrix(M, nrows, ncols)
        implicit none

        integer*4 ii, nrows, ncols
        real*8 M(nrows, ncols)

        do ii = 1, nrows
            print '(20f7.2)', M(ii, 1:ncols)
        end do
    end subroutine

    ! print square complex matrix in a nice format
    !
    ! Inputs:
    !   M: Matrix to print
    !   n: Number of rows and columns
    !
    subroutine print_complex_matrix(M, n)
        implicit none

        integer n
        integer ii
        complex*16 M(n, n)

        do ii = 1, n
            print '(*(sp, f7.4, 1x, f7.4, "i", 3x))', M(ii, :)
        end do
    end subroutine

end module
