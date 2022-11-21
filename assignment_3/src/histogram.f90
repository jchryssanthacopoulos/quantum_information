! this module contains a subroutine writes a histogram to a file
module histogram

contains
    ! creates a histogram from set of values and writes it to a file
    !
    ! Inputs:
    !   values: Values to bin into a histogram
    !   nbins: Number of bins to use
    !   min_val: Minimum value of the histogram
    !   max_val: Maximum value of the histogram
    !   filename: Name of file to save
    !
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
        write(1, *) "left_edges,centers,count,norm_count"

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
            write(1, '(f10.7, ",", f10.7, ",", i10, ",", f10.7)') &
                left_bin_edge, left_bin_edge + dx / 2, bin_count, bin_count / norm_factor

            left_bin_edge = left_bin_edge + dx
        end do

        close(1)
    end subroutine

end module
