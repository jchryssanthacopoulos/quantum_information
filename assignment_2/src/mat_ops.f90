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

end module
