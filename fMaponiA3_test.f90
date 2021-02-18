program Interface_test
    use Sherman_Morrison, only : MaponiA3
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    implicit none

    integer i, j, col !! Iterators
    integer(c_int) :: Dim, N_updates
    integer(c_int), dimension(:), allocatable :: Updates_index
    real(c_double), dimension(:,:), allocatable :: S, A, Updates
    real(c_double), dimension(:,:), allocatable :: S_inv, A_inv

    Dim = 4
    N_updates = 2
    allocate( S(Dim,Dim), S_inv(Dim,Dim), A(Dim,Dim), A_inv(Dim,Dim), &
              Updates(Dim,N_updates), Updates_index(N_updates))

    !! Initialize S, S_inv, A and A_inv
    S(1,1) = 1.0d0
    S(1,2) = 0.0d0
    S(1,3) = 1.0d0
    S(1,4) = -1.0d0
    S(2,1) = 0.0d0
    S(2,2) = 1.0d0
    S(2,3) = 1.0d0
    S(2,4) = 0.0d0
    S(3,1) = -1.0d0
    S(3,2) = 0.0d0
    S(3,3) = -1.0d0
    S(3,4) = 0.0d0
    S(4,1) = 1.0d0
    S(4,2) = 1.0d0
    S(4,3) = 1.0d0
    S(4,4) = 1.0d0

    S_inv(1,1) = 1.0d0
    S_inv(1,2) = -1.0d0
    S_inv(1,3) = 1.0d0
    S_inv(1,4) = 1.0d0
    S_inv(2,1) = 1.0d0
    S_inv(2,2) = 0.0d0
    S_inv(2,3) = 2.0d0
    S_inv(2,4) = 1.0d0
    S_inv(3,1) = -1.0d0
    S_inv(3,2) = 1.0d0
    S_inv(3,3) = -2.0d0
    S_inv(3,4) = -1.0d0
    S_inv(4,1) = -1.0d0
    S_inv(4,2) = 0.0d0
    S_inv(4,3) = -1.0d0
    S_inv(4,4) = 0.0d0

    A(1,1) = 1.0d0
    A(1,2) = 0.0d0
    A(1,3) = 1.0d0
    A(1,4) = -1.0d0
    A(2,1) = 0.0d0
    A(2,2) = -1.0d0
    A(2,3) = 1.0d0
    A(2,4) = -1.0d0
    A(3,1) = -1.0d0
    A(3,2) = 0.0d0
    A(3,3) = -1.0d0
    A(3,4) = 0.0d0
    A(4,1) = 1.0d0
    A(4,2) = 1.0d0
    A(4,3) = 1.0d0
    A(4,4) = 1.0d0

    A_inv(1,1) = 0.0d0
    A_inv(1,2) = -1.0d0
    A_inv(1,3) = -2.0d0
    A_inv(1,4) = -1.0d0
    A_inv(2,1) = 1.0d0
    A_inv(2,2) = 0.0d0
    A_inv(2,3) = 2.0d0
    A_inv(2,4) = 1.0d0
    A_inv(3,1) = 0.0d0
    A_inv(3,2) = 1.0d0
    A_inv(3,3) = 1.0d0
    A_inv(3,4) = 1.0d0
    A_inv(4,1) = -1.0d0
    A_inv(4,2) = 0.0d0
    A_inv(4,3) = -1.0d0
    A_inv(4,4) = 0.0d0

    !! Prepare Updates and Updates_index
    Updates(1,1) = 0
    Updates(1,2) = 0
    Updates(2,1) = -2
    Updates(2,2) = -1
    Updates(3,1) = 0
    Updates(3,2) = 0
    Updates(4,1) = 0
    Updates(4,2) = 0

    Updates_index(1) = 2
    Updates_index(2) = 4

    write(*,*)
    write(*,*) "Old S = "
    do i=1,Dim
        do j=1,Dim
            write(*,"(F3.0,3X)", advance="no") S(i,j)
        end do
        write(*,*)
    end do

    write(*,*)
    write(*,*) "Old S_inv = "
    do i=1,Dim
        do j=1,Dim
            write(*,"(F3.0,3X)", advance="no") S_inv(i,j)
        end do
        write(*,*)
    end do

    write(*,*)
    write(*,*) "Updates = "
    do i=1,Dim
        do j=1,N_updates
            write(*,"(F3.0,3X)", advance="no") Updates(i,j)
        end do
        write(*,*)
    end do

    write(*,*)
    write(*,*) "Updates_index = "
    do i=1,N_updates
            write(*,"(I1,3X)", advance="no") Updates_index(i)
    end do
    write(*,*)

    !! Update S
    do i=1,N_updates
        do j=1,Dim
            col = Updates_index(i)
            S(j,col) = S(j,col) + Updates(j,i)
        end do
    end do
            
    !! Update S_inv
    call MaponiA3(S_inv, Dim, N_updates, Updates, Updates_index)

    write(*,*)
    write(*,*)
    write(*,*) "New computed S = "
    do i=1,Dim
        do j=1,Dim
            write(*,"(F3.0,3X)", advance="no") S(i,j)
        end do
        write(*,*)
    end do

    write(*,*)
    write(*,*) "New actual S = "
    do i=1,Dim
        do j=1,Dim
            write(*,"(F3.0,3X)", advance="no") A(i,j)
        end do
        write(*,*)
    end do

    write(*,*)
    write(*,*)
    write(*,*) "New computed S_inv = "
    do i=1,Dim
        do j=1,Dim
            write(*,"(F3.0,3X)", advance="no") S_inv(i,j)
        end do
        write(*,*)
    end do
    
    write(*,*)
    write(*,*) "New actual S_inv = "
    do i=1,Dim
        do j=1,Dim
            write(*,"(F3.0,3X)", advance="no") A_inv(i,j)
        end do
        write(*,*)
    end do

    deallocate(Updates_index, A, A_inv, S, Updates, S_inv)
end program
