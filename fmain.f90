program Interface_test
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    use MYMODULE, only : MYSUBROUTINE
    implicit none

    integer i, j !! Iterators
    integer(c_int) :: dim, N_updates
    integer(c_int), dimension(:), allocatable :: Ar_index
    integer(c_int), dimension(:,:), allocatable :: A, A0, Ar
    real(c_double), dimension(:,:), allocatable :: A0_inv

    dim = 3
    N_updates = dim
    allocate(Ar_index(dim), A(dim,dim), A0(dim,dim), Ar(dim,dim), A0_inv(dim,dim))

    !! Initialize A with M=3 and fill acc. to Eq. (17) from paper
    A(1,1) = 1
    A(1,2) = 1
    A(1,3) = -1
    A(2,1) = 1
    A(2,2) = 1
    A(2,3) = 0
    A(3,1) = -1
    A(3,2) = 0
    A(3,3) = -1

    do i=1,3
        do j=1,3
            write(*,"(I)", advance="no") A(i,j)
        end do
        write(*,*)
    end do

    deallocate(Ar_index, A, A0, Ar, A0_inv)
end program
