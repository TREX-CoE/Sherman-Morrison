program QMCChem_dataset_test
    use Sherman_Morrison
    use Helpers
    implicit none

    integer :: i, j, col !! Iterators
    integer :: cycle_id, dim, n_updates
    integer(c_int), dimension(:), allocatable :: Updates_index
    real(c_double), dimension(:,:), allocatable :: Updates
    real(c_double), dimension(:,:), allocatable :: S, S_inv, S_inv_t

    call Read_dataset("datasets/update_cycle_13.dat", &
                       cycle_id, &
                       dim, &
                       n_updates, &
                       S, &
                       S_inv, &
                       Updates_index, &
                       Updates)
    allocate(S_inv_t(dim,dim))

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
    do j=1,n_updates
        do i=1,dim
            col = Updates_index(j)
            S(i,col) = S(i,col) + Updates(i,j)
        end do
    end do

    !! Update S_inv
    !! S_inv needs to be transposed first before it
    !! goes to MaponiA3
    call Transpose(S_inv, S_inv_t, dim)
    call MaponiA3(S_inv_t, dim, n_updates, Updates, Updates_index)
    !! S_inv_t needs to be transposed back before we
    !! can multiply it with S to test unity
    call Transpose(S_inv_t, S_inv, dim)

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

    deallocate(S, S_inv, S_inv_t, Updates, Updates_index)
end program