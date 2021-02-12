program QMCChem_dataset_test
    use Sherman_Morrison, only : MaponiA3
    use Utils, only : Read_dataset
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    implicit none

    integer :: i, j, col !! Iterators
    integer :: cycle_id, dim, n_updates
    integer(c_int), dimension(:), allocatable :: Updates_index
    real(c_double), dimension(:,:), allocatable :: S, S_inv, Updates

    call Read_dataset("test.dataset.dat", &
                       cycle_id, &
                       dim, &
                       n_updates, &
                       S, &
                       S_inv, &
                       Updates_index, &
                       Updates)

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
