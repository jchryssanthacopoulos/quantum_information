program tutorial_2
    implicit none

    print *, "Square of 5 is", mysquare(5d0)
    print *, "5 - 2 is", my_subtract_2(5d0)

contains
    ! define functions at end
    ! example using result
    function mysquare(x) result(xs)
        implicit none
        real*8 x
        real*8 xs
        xs = x ** 2
    end function

    ! example using return
    function my_subtract_2(x)
        implicit none
        real*8 x
        real*8 my_subtract_2
        my_subtract_2 = x - 2
        return
    end function

end program tutorial_2
