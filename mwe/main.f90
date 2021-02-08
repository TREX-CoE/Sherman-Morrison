program test
  use iso_c_binding
  implicit none

  interface
    subroutine hello() bind(C, name="worker")
    end subroutine
  end interface

  call hello()

end program test
