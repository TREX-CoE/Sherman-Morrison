program Interface_test
    use Sherman_Morrison, only : MaponiA3
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    implicit none

    integer i, j !! Iterators
    integer(c_int) :: Dim, N_updates
    integer(c_int), dimension(:), allocatable :: Updates_index
    real(c_double), dimension(:,:), allocatable :: A, S, Updates
    real(c_double), dimension(:,:), allocatable :: S_inv

    Dim = 3
    N_updates = 3
    allocate(Updates_index(Dim), A(Dim,Dim), S(Dim,Dim), Updates(Dim,Dim), S_inv(Dim,Dim))

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

    !! Prepare the diagonal matrix S and the update matrix Updates
    do i=1,Dim
        Updates_index(i) = i
        do j=1,Dim
            if (i == j) then
                S(i,j) = A(i,j)
                S_inv(i,j) = 1.0d0 / S(i,j)
            else
                S(i,j) = 0.0d0
                S_inv(i,j) = 0.0d0
            end if
            Updates(i,j) = A(i,j) - S(i,j)
        end do
    end do

    call MaponiA3(S_inv, Dim, N_updates, Updates, Updates_index)

    do i=1,Dim
        do j=1,Dim
            write(*,"(F3.0,3X)", advance="no") S_inv(i,j)
        end do
        write(*,*)
    end do

    deallocate(Updates_index, A, S, Updates, S_inv)
end program
