module Utils
    implicit none

    contains
        subroutine Read_dataset(filename, cycle_id, dim, n_updates &
                    S, S_inv, Updates_index, Updates)
            
            character (len = 64), intent(in) :: filename
            
            
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
        end subroutine Read_dataset
    end module Utils


end module Sherman_Morrison
