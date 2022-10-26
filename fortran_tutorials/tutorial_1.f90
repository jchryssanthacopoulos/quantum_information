! define modules
module intro
    integer*2 ii
    real*8 test1
end module intro


module advanced
    ! real*8 = double precision
    ! :: separates declaration from name
    double precision, dimension(:), allocatable :: vect

    ! include 'library.f'

    contains
        function mysqrt(x) result(sx)
            real*8 x
            real*8 sx
            sx = sqrt(x)
        end function
end module advanced


program test
    use intro  ! import modules
    use advanced
    implicit none  ! means there are no implicit declarations

    real*8 sx
    integer y

    ! allocate memory for the vector
    allocate(vect(10))

    ! for loop
    do ii = 1, 10
        vect(ii) = ii
    end do

    sx = mysqrt(5d0)

    ! * represents format (default)
    print *, "Square root of 5 is", sx
    print *, "5 squared is", mysquare(5d0)

    ! call myopen("my_file.txt")

    ! deallocate
    deallocate(vect)

    y = 10
    if (y.gt.0) then
        print *, "y is greater than 0!"
    else
        print *, "y is not greater than 0"
    endif

    stop

contains
    function mysquare(x) result(xs)
        real*8 x
        real*8 xs
        xs = x ** 2
    end function

end program


subroutine myopen(file1)
    ! use advanced
    implicit none

    character*11 file1
    file1 = trim(file1)
    print *, "File = ", file1

    open(unit = 50, file = file1, status = 'unknown')
    write(50, *) "Hello!"
    close(50)

    return
end subroutine
