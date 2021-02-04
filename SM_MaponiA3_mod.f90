module MYMODULE
    interface
        subroutine MYSUBROUTINE(Slater0, Slater_inv, dim, n_updates, Updates, Updates_index) bind(C, name="Sherman_Morrison")
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in) :: dim, n_updates
            integer(c_int), dimension(:), allocatable, intent(in) :: Updates_index
            integer(c_int), dimension(:,:), allocatable, intent(in) :: Slater0, Updates
            real(c_double), dimension(:,:), allocatable, intent(in out) :: Slater_inv
        end subroutine MYSUBROUTINE
    end interface
end module MYMODULE