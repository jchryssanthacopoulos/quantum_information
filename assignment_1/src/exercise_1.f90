! define modules
module first_module
    implicit none

    real*8 var1, var2

    contains
        function mysqrt(x) result(sx)
            real*8 x
            real*8 sx
            sx = sqrt(x)
        end function

end module first_module


program exercise_1
    use first_module  ! import modules
    implicit none  ! means there are no implicit declarations

    var1 = 5d0
    var2 = mysqrt(var1)

    print "('The square root of ', f3.1, ' is ', f6.4)", var1, var2

end program
