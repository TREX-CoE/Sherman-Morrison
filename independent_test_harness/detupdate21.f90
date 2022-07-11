subroutine detupd(vDim, vLDS, Updates, Updates_index, &
  Slater_inv, determinant) bind(C)
  use iso_c_binding
  implicit none                                                         ! det.irp.f_template_577: 428

  external :: det_update21
  integer(c_int64_t), intent(in), value :: vLDS, vDim
  real(c_double), intent(in)   :: Updates(vLDS)
  integer(c_int64_t), intent(in)  :: Updates_index(1)
  real(c_double), intent(inout) :: Slater_inv(vLDS,vDim)
  real(c_double), intent(inout) :: determinant

  integer(c_int64_t) :: l, n, LDS                                   ! det.irp.f_template_577: 432

  n = vDim
  LDS = vLDS
  l = Updates_index(1)
  call det_update21(n, LDS, Updates, l, Slater_inv, determinant)      ! det.irp.f_template_577: 427
end

subroutine det_update21(n, LDS, u, l, S_inv, d)                 ! det.irp.f_template_577: 427
  use iso_c_binding
  implicit none                                                         ! det.irp.f_template_577: 428
  integer(c_int64_t), intent(in) :: l                               ! det.irp.f_template_577: 430
  integer(c_int64_t), intent(in) :: n,LDS                               ! det.irp.f_template_577: 430
  real(c_double),intent(inout) :: S_inv(LDS,n)                       ! det.irp.f_template_577: 435
  real(c_double),intent(inout) :: d                                   ! det.irp.f_template_577: 436
  real(c_double), intent(in)   :: u(n)
  real(c_double)               :: z(n), w(n), lambda, d_inv  ! det.irp.f_template_577: 438
  integer(c_int64_t)           :: i,j                                 ! det.irp.f_template_577: 444
  real(c_double)               :: zj, zj1, zj2, zj3                  ! det.irp.f_template_577: 445

  !DIR$ ATTRIBUTES ALIGN : 32 :: z, w                                ! det.irp.f_template_577: 439
  !DIR$ ASSUME_ALIGNED u : 32                                           ! det.irp.f_template_577: 440
  !DIR$ ASSUME_ALIGNED S_inv : 32                                       ! det.irp.f_template_577: 441
  !DIR$ ASSUME (mod(LDS,32/8) == 0)                                     ! det.irp.f_template_577: 442
  !DIR$ ASSUME (LDS >= 21)                                              ! det.irp.f_template_577: 443
  
  zj = 0.d0    !! dot prod col S_inv and update: vT*S_inv*u                                                       ! det.irp.f_template_577: 451
  !DIR$ NOPREFETCH                                                      ! det.irp.f_template_577: 452
  do i=1,21-1,4                                                         ! det.irp.f_template_577: 453
    zj = zj + S_inv(i,l)*u(i) + S_inv(i+1,l)*u(i+1)  &
            + S_inv(i+2,l)*u(i+2) + S_inv(i+3,l)*u(i+3) ! det.irp.f_template_577: 454
  enddo                                                                 ! det.irp.f_template_577: 456
  zj = zj + S_inv(21,l)*u(21)                                           ! det.irp.f_template_577: 457

  d_inv = 1.d0/d  ! reciprocal of old det                                                      ! det.irp.f_template_577: 459
  d = d+zj                                                              ! det.irp.f_template_577: 460
  lambda = d*d_inv                                                      ! det.irp.f_template_577: 461
  if (dabs(lambda) < 1.d-3) then                                        ! det.irp.f_template_577: 462
    ! d = 0.d0                                                            ! det.irp.f_template_577: 463
    return                                                              ! det.irp.f_template_577: 464
  endif                                                                 ! det.irp.f_template_577: 465
  !DIR$ VECTOR ALIGNED                                                  ! det.irp.f_template_577: 467
  do j=1,21-1,4                                                         ! det.irp.f_template_577: 468
   zj  = 0.d0                                                           ! det.irp.f_template_577: 469
   zj1 = 0.d0                                                           ! det.irp.f_template_577: 470
   zj2 = 0.d0                                                           ! det.irp.f_template_577: 471
   zj3 = 0.d0                                                           ! det.irp.f_template_577: 472
   !DIR$ VECTOR ALIGNED                                                 ! det.irp.f_template_577: 473
   !DIR$ NOPREFETCH                                                     ! det.irp.f_template_577: 474
   do i=1,21-1                                                          ! det.irp.f_template_577: 475
    zj  = zj  + S_inv(i,j  )*u(i)                                       ! det.irp.f_template_577: 476
    zj1 = zj1 + S_inv(i,j+1)*u(i)                                       ! det.irp.f_template_577: 477
    zj2 = zj2 + S_inv(i,j+2)*u(i)                                       ! det.irp.f_template_577: 478
    zj3 = zj3 + S_inv(i,j+3)*u(i)                                       ! det.irp.f_template_577: 479
   enddo                                                                ! det.irp.f_template_577: 480
   z(j  ) = zj  + S_inv(21,j  )*u(21)                                   ! det.irp.f_template_577: 481
   z(j+1) = zj1 + S_inv(21,j+1)*u(21)                                   ! det.irp.f_template_577: 482
   z(j+2) = zj2 + S_inv(21,j+2)*u(21)                                   ! det.irp.f_template_577: 483
   z(j+3) = zj3 + S_inv(21,j+3)*u(21)                                   ! det.irp.f_template_577: 484
  enddo                                                                 ! det.irp.f_template_577: 485
  zj  = 0.d0                                                            ! det.irp.f_template_577: 487
  !DIR$ VECTOR ALIGNED                                                  ! det.irp.f_template_577: 488
  !DIR$ NOPREFETCH                                                      ! det.irp.f_template_577: 489
  do i=1,21-1                                                           ! det.irp.f_template_577: 490
   zj = zj + S_inv(i,21)*u(i)                                           ! det.irp.f_template_577: 491
  enddo                                                                 ! det.irp.f_template_577: 492
  z(21) = zj + S_inv(21,21)*u(21)                                       ! det.irp.f_template_577: 493
  !DIR$ NOPREFETCH                                                      ! det.irp.f_template_577: 495
  do i=1,21                                                             ! det.irp.f_template_577: 496
    w(i) = S_inv(i,l)*d_inv                                             ! det.irp.f_template_577: 497
  enddo                                                                 ! det.irp.f_template_577: 499
  do i=1,21-1,4                                                         ! det.irp.f_template_577: 501
   zj  = z(i  )                                                         ! det.irp.f_template_577: 502
   zj1 = z(i+1)                                                         ! det.irp.f_template_577: 503
   zj2 = z(i+2)                                                         ! det.irp.f_template_577: 504
   zj3 = z(i+3)                                                         ! det.irp.f_template_577: 505
   !DIR$ VECTOR ALIGNED                                                 ! det.irp.f_template_577: 506
   !DIR$ NOPREFETCH                                                     ! det.irp.f_template_577: 507
   do j=1,21-1                                                          ! det.irp.f_template_577: 508
    S_inv(j,i  ) = S_inv(j,i  )*lambda - w(j)*zj                        ! det.irp.f_template_577: 509
    S_inv(j,i+1) = S_inv(j,i+1)*lambda - w(j)*zj1                       ! det.irp.f_template_577: 510
    S_inv(j,i+2) = S_inv(j,i+2)*lambda - w(j)*zj2                       ! det.irp.f_template_577: 511
    S_inv(j,i+3) = S_inv(j,i+3)*lambda - w(j)*zj3                       ! det.irp.f_template_577: 512
   enddo                                                                ! det.irp.f_template_577: 513
   S_inv(21,i  ) = S_inv(21,i  )*lambda - w(21)*zj                      ! det.irp.f_template_577: 514
   S_inv(21,i+1) = S_inv(21,i+1)*lambda - w(21)*zj1                     ! det.irp.f_template_577: 515
   S_inv(21,i+2) = S_inv(21,i+2)*lambda - w(21)*zj2                     ! det.irp.f_template_577: 516
   S_inv(21,i+3) = S_inv(21,i+3)*lambda - w(21)*zj3                     ! det.irp.f_template_577: 517
  enddo                                                                 ! det.irp.f_template_577: 518
  zj = z(21)                                                            ! det.irp.f_template_577: 520
  !DIR$ VECTOR ALIGNED                                                  ! det.irp.f_template_577: 521
  !DIR$ NOPREFETCH                                                      ! det.irp.f_template_577: 522
  do i=1,21                                                             ! det.irp.f_template_577: 523
   S_inv(i,21) = S_inv(i,21)*lambda -w(i)*zj                            ! det.irp.f_template_577: 524
  enddo                                                                 ! det.irp.f_template_577: 525
end                                                                     ! det.irp.f_template_577: 528
