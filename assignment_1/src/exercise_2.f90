program exercise_2
    implicit none

    integer*2 int1_2, int2_2
    integer*4 int1_4, int2_4
    real*4 real1_4, real2_4
    real*8 real1_8, real2_8

    ! sum using integer*2
    int1_2 = 2000000
    int2_2 = 1
    print "('The sum of ', i6, ' and ', i1, ' using INTEGER*2 is ', i6)", int1_2, int2_2, int1_2 + int2_2

    ! sum using integer*4
    int1_4 = 2000000
    int2_4 = 1
    print "('The sum of ', i7, ' and ', i1, ' using INTEGER*4 is ', i7)", int1_4, int2_4, int1_4 + int2_4

    ! sum using real*4
    real1_4 = 4e0*atan(1e0)*1e32
    real2_4 = sqrt(2e0)*1e21
    print "('The sum of ', es14.8, ' and ', es14.8, ' using REAL*4 is ', es14.8)", real1_4, real2_4, real1_4 + real2_4

    ! sum using real*8
    real1_8 = 4d0*datan(1d0)*1d32
    real2_8 = sqrt(2d0)*1d21
    print "('The sum of ', es22.16, ' and ', es22.16, ' using REAL*8 is ', es22.16)", real1_8, real2_8, real1_8 + real2_8

end program
