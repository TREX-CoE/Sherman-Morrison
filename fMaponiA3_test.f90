program Interface_test
    use Sherman_Morrison, only : MaponiA3
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    implicit none

    integer i, j !! Iterators
    integer(c_int) :: dim, N_updates
    integer(c_int), dimension(:), allocatable :: Ar_index
    real(c_double), dimension(:,:), allocatable :: A, A0, Ar
    real(c_double), dimension(:,:), allocatable :: A0_inv

    dim = 3
    N_updates = 3
    allocate(Ar_index(dim), A(dim,dim), A0(dim,dim), Ar(dim,dim), A0_inv(dim,dim))

    !! Initialize A with M=3 and fill acc. to Eq. (17) from paper
    A(1,1) = 1.0d0
    A(1,2) = 1.0d0
    A(1,3) = -1.0d0
    A(2,1) = 1.0d0
    A(2,2) = 1.0d0
    A(2,3) = 0.0d0
    A(3,1) = -1.0d0
    A(3,2) = 0.0d0
    A(3,3) = -1.0d0

    !! Prepare the diagonal matrix A0 and the update matrix Ar
    do i=1,dim
        Ar_index(i) = i
        do j=1,dim
            if (i == j) then
                A0(i,j) = A(i,j)
                A0_inv(i,j) = 1.0d0 / A0(i,j)
            else
                A0(i,j) = 0.0d0
                A0_inv(i,j) = 0.0d0
            end if
            Ar(i,j) = A(i,j) - A0(i,j)
        end do
    end do

    call MaponiA3(A0, A0_inv, dim, n_updates, Ar, Ar_index)

    do i=1,dim
        do j=1,dim
            write(*,"(F3.0,3X)", advance="no") A0_inv(i,j)
        end do
        write(*,*)
    end do

    deallocate(Ar_index, A, A0, Ar, A0_inv)
end program
