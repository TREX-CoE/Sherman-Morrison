program Interface_test
    use Sherman_Morrison
    use Helpers
    implicit none

    integer i, j, col !! Iterators
    integer(c_int) :: Dim, N_updates
    integer(c_int), dimension(:), allocatable :: Updates_index
    real(c_double), dimension(:,:), allocatable :: S, A, Updates
    real(c_double), dimension(:,:), allocatable :: S_inv, S_inv_t, A_inv

    Dim = 4
    N_updates = 2
    allocate( S(Dim,Dim), S_inv(Dim,Dim), S_inv_t(Dim,Dim), A(Dim,Dim), A_inv(Dim,Dim), &
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

    !! Write current S and S_inv to file for check in Octave
    open(unit = 2000, file = "Slater_old.dat")
    open(unit = 3000, file = "Slater_old_inv.dat")
    do i=1,dim
        do j=1,dim
            write(2000,"(E23.15, 1X)", advance="no") S(i,j)
            write(3000,"(E23.15, 1X)", advance="no") S_inv(i,j)
        end do
        write(2000,*)
        write(3000,*)
    end do
    close(2000)
    close(3000)

    !! Write Updates to file to check
    open(unit = 2000, file = "Updates.dat")
    do i=1,dim
        do j=1,n_updates
            write(2000,"(E23.15, 1X)", advance="no") Updates(i,j)
        end do
        write(2000,*)
    end do
    close(2000)

    !! Update S
    do i=1,N_updates
        do j=1,Dim
            col = Updates_index(i)
            S(j,col) = S(j,col) + Updates(j,i)
        end do
    end do
            
    !! Update S_inv
    call Transpose(S_inv, S_inv_t, Dim)
    call MaponiA3(S_inv_t, Dim, N_updates, Updates, Updates_index)
    call Transpose(S_inv_t, S_inv, Dim)

    !! Write new S and S_inv to file for check in Octave
    open(unit = 4000, file = "Slater.dat")
    open(unit = 5000, file = "Slater_inv.dat")
    do i=1,dim
        do j=1,dim
            write(4000,"(E23.15, 1X)", advance="no") S(i,j)
            write(5000,"(E23.15, 1X)", advance="no") S_inv(i,j)
        end do
        write(4000,*)
        write(5000,*)
    end do
    close(4000)
    close(5000)

    deallocate(Updates_index, A, A_inv, S, Updates, S_inv)
end program
