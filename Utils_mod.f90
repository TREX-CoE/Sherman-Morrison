module Utils
    implicit none
    contains
        subroutine Read_dataset(filename, &
                                cycle_id, &
                                dim, &
                                n_updates, &
                                S, &
                                S_inv, &
                                Updates_index, &
                                Updates)
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            implicit none
            
            character (len = *), intent(in) :: filename
            integer, intent(inout) :: cycle_id, dim, n_updates
            integer(c_int), allocatable, intent(inout) :: Updates_index(:)
            real(c_double), allocatable, intent(inout) :: S(:,:), S_inv(:,:)
            real(c_double), allocatable, intent(inout) :: Updates(:,:)
            integer :: i, j
            character (len = 32) :: ignore

            !! Start reading dataset from file
            open(unit = 1000, file = filename)
            read(1000,*)
            read(1000,*) ignore, cycle_id
            read(1000,*) ignore, dim
            read(1000,*) ignore,n_updates

            !! Initialise the arrays
            allocate(Updates_index(n_updates), &
                     S(dim,dim), &
                     S_inv(dim,dim), &
                     Updates(dim,n_updates))
            
            !! Read S and S_inv
            read(1000,*)
            do i=1,dim
                do j=1,dim
                    !! For some reason S_inv needs to be read transposed &
                    !  for S*S_inv to give the identity matrix
                    read(1000,*) ignore, ignore, S(i,j), S_inv(j,i)
                end do
            end do

            !! Read Updates and Updates_index
            do j=1,n_updates
                read(1000,*) ignore, Updates_index(j)
                do i=1,dim
                    read(1000,*) ignore, Updates(i,j)
                end do
            end do
            !! End of reading the dataset from file
           close(1000)
        end subroutine Read_dataset
end module Utils