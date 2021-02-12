program QMCChem_dataset_test
    use Sherman_Morrison, only : MaponiA3
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    implicit none

    integer :: i, j, col !! Iterators
    integer :: cycle_id, dim, n_updates
    integer(c_int), dimension(:), allocatable :: Updates_index
    real(c_double), dimension(:,:), allocatable :: S, S_inv, Updates
    character (len = 32) :: ignore

    !! Start of reading the dataset from file
    open(unit = 1000, file = "test.dataset.dat")
    read(1000,*)
    read(1000,*) ignore, cycle_id
    read(1000,*) ignore, dim
    read(1000,*) ignore,n_updates
   
    allocate(Updates_index(n_updates), S(dim,dim), &
                S_inv(dim,dim), Updates(dim,n_updates))
    
    !! Read S and S_inv
    read(1000,*)
    do i=1,dim
        do j=1,dim
            read(1000,*) ignore, ignore, S(i,j), S_inv(i,j)
        end do
    end do

    !! Read the updates Updates and Updates_index
    do j=1,n_updates
        read(1000,*) ignore, Updates_index(j)
        do i=1,dim
            read(1000,*) ignore, Updates(i,j)
        end do
    end do

    read(1000,*) ignore
    close(1000)
    !! End of reading the dataset from file

    !! Write current S and S_inv to file for check in Octave
    open(unit = 2000, file = "Slater_old.dat")
    open(unit = 3000, file = "Slater_inv_old.dat")
    do i=1,dim
        do j=1,dim
            write(2000,"(E23.15, 1X)", advance="no") S(i,j)
            write(3000,"(E23.15, 1X)", advance="no") S_inv(i,j)
        end do
        write(2000,*)
        write(3000,*)
    end do
    call flush(2000)
    call flush(3000)
    close(2000)
    close(3000)

    !! Send S, S_inv and Updates to MaponiA3 algo
    call MaponiA3(S, S_inv, dim, n_updates, Updates, Updates_index)

    !! Update S itself
    do j=1,n_updates
        do i=1,dim
            col = Updates_index(j)
            S(i,col) = S(i,col) + Updates(i,col)
        end do
    end do

    !! Write new S and S_inv to file for check in Octave
    open(unit = 2000, file = "Slater_new.dat")
    open(unit = 3000, file = "Slater_inv_new.dat")
    do i=1,dim
        do j=1,dim
            write(2000,"(E23.15, 1X)", advance="no") S(i,j)
            write(3000,"(E23.15, 1X)", advance="no") S_inv(i,j)
        end do
        write(2000,*)
        write(3000,*)
    end do
    call flush(2000)
    call flush(3000)

    close(2000)
    close(3000)

    deallocate(S, S_inv, Updates, Updates_index)
end program
