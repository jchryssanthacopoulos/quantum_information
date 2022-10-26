MODULE genswap
    IMPLICIT NONE

    TYPE point
        REAL :: x, y
    END TYPE point

    INTERFACE swap ! generic interface
    MODULE PROCEDURE swappoint
    END INTERFACE

    CONTAINS
        SUBROUTINE swappoint (a, b)
            TYPE (point), INTENT(INOUT) :: a, b
            TYPE (point) :: temp
            temp = a; a = b; b = temp
        END SUBROUTINE swappoint

END MODULE genswap


! module matrices
!     type dmatrix
!         ! define components of the type
!         integer, dimension(2) :: N  ! vector with two integers
!         double complex, dimension(:, :), allocatable :: elem
!     end type dmatrix

!     interface operator (.Adj.)
!     module procedure MatAdjoint, VecAdjoint
!     end interface

!     contains
!         function MatAdjoint(x)
!             print *, "A"
!         end function
! end module matrices


program tutorial_2
    ! use matrices
    use genswap
    implicit none

    type (point) :: a, b
    ! type (dmatrix) :: a, b

    ! give scalar to all elements of matrix
    ! a%elem = (0d0, 1d0)
    ! b = .Adj.a

    a%x = 2
    a%y = 3

    b%x = 4
    b%y = 5

    call swappoint(a, b)

    print *, a%x, a%y
    print *, b%x, b%y

end program tutorial_2
