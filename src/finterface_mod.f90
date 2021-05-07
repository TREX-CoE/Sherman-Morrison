module Sherman_Morrison
    use, intrinsic :: iso_c_binding, only : c_int, c_double
    interface
        subroutine MaponiA3(Slater_inv, dim, n_updates, Updates, Updates_index) bind(C, name="MaponiA3_f")
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in) :: dim, n_updates
            integer(c_int), dimension(:), allocatable, intent(in) :: Updates_index
            real(c_double), dimension(:,:), allocatable, intent(in) :: Updates
            real(c_double), dimension(:,:), allocatable, intent(in out) :: Slater_inv
        end subroutine MaponiA3
        subroutine MaponiA3S(Slater_inv, dim, n_updates, Updates, Updates_index) bind(C, name="MaponiA3S_f")
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in) :: dim, n_updates
            integer(c_int), dimension(:), allocatable, intent(in) :: Updates_index
            real(c_double), dimension(:,:), allocatable, intent(in) :: Updates
            real(c_double), dimension(:,:), allocatable, intent(in out) :: Slater_inv
        end subroutine MaponiA3S
        subroutine SM1(Slater_inv, dim, n_updates, Updates, Updates_index) bind(C, name="SM1_f")
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in) :: dim, n_updates
            integer(c_int), dimension(:), allocatable, intent(in) :: Updates_index
            real(c_double), dimension(:,:), allocatable, intent(in) :: Updates
            real(c_double), dimension(:,:), allocatable, intent(in out) :: Slater_inv
        end subroutine SM1
        subroutine SM2(Slater_inv, dim, n_updates, Updates, Updates_index) bind(C, name="SM2_f")
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in) :: dim, n_updates
            integer(c_int), dimension(:), allocatable, intent(in) :: Updates_index
            real(c_double), dimension(:,:), allocatable, intent(in) :: Updates
            real(c_double), dimension(:,:), allocatable, intent(in out) :: Slater_inv
        end subroutine SM2
        subroutine SM3(Slater_inv, dim, n_updates, Updates, Updates_index) bind(C, name="SM3_f")
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in) :: dim, n_updates
            integer(c_int), dimension(:), allocatable, intent(in) :: Updates_index
            real(c_double), dimension(:,:), allocatable, intent(in) :: Updates
            real(c_double), dimension(:,:), allocatable, intent(in out) :: Slater_inv
        end subroutine SM3
        subroutine SM4(Slater_inv, dim, n_updates, Updates, Updates_index) bind(C, name="SM4_f")
            use, intrinsic :: iso_c_binding, only : c_int, c_double
            integer(c_int), intent(in) :: dim, n_updates
            integer(c_int), dimension(:), allocatable, intent(in) :: Updates_index
            real(c_double), dimension(:,:), allocatable, intent(in) :: Updates
            real(c_double), dimension(:,:), allocatable, intent(in out) :: Slater_inv
        end subroutine SM4
    end interface
end module Sherman_Morrison
