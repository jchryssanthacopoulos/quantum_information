! Exercise 3
! ==========
!
! This program creates a random matrix of double complex numbers with dimensions provided by the user, computes
!   the trace and adjoint, and writes the matrices to files
!
! Inputs:
!   number of rows of matrix
!   number of columns of matrix
!
! Returns:
!   The created matrix, its trace and adjoint
!   It also writes the original matrix and adjoint to files mat.txt and mat_adj.txt, respectively
!


program exercise_3
    use arg_parse
    use cmatrix_type
    implicit none

    character(len=20) char_input(2)
    integer*4 dims(2)
    type(cmatrix) M, Mdagger

    print *, "Enter the number of rows and columns of your matrix:"
    read *, char_input

    ! check if dimensions are integers
    call check_dims_integers(char_input, dims, 2)
    if (status /= 0) then
        print *, "Dimensions need to be integers!"
        stop
    end if

    ! check if dimensions are positive
    call check_dims_positive(dims, 2)
    if (status /= 0) then
        print *, "Dimensions must be greater than zero!"
        stop
    end if

    ! initialize matrix with random data
    M = Init(dims(1), dims(2))

    ! display matrix and save to file
    print *, "The original matrix is:"
    call print_matrix(M)
    call print_trace(M)

    print *, "Saving to file mat.txt ..."
    call write_matrix(M, "mat.txt")

    ! get and display adjoint and save to file
    Mdagger = .Adj.M
    print *, "The adjoint matrix is:"
    call print_matrix(Mdagger)
    call print_trace(Mdagger)

    print *, "Saving to file mat_adj.txt ..."
    call write_matrix(Mdagger, "mat_adj.txt")

    print *, "Deleting matrices ..."
    call Del(M)
    call Del(Mdagger)

end program
