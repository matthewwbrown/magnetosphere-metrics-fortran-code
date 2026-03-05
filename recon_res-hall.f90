program main

  use inputOutput, only : save_rate_to_netcdf, load_field_from_netcdf

  implicit none

  real(8), dimension(:,:,:), allocatable :: Bx,By,Bz,magB
  real(8), dimension(:,:,:), allocatable :: Fx,Fy,Fz,magF
  real(8), dimension(:,:,:), allocatable :: Jx,Jy,Jz
  real(8), dimension(:,:,:), allocatable :: Bfx,Bfy,Bfz
  real(8), dimension(:,:,:), allocatable :: lambda,alpha
  real(8), dimension(:,:,:), allocatable :: rate

  ! Hall-related
  real(8), dimension(:,:,:), allocatable :: JxB_x,JxB_y,JxB_z
  real(8), dimension(:,:,:), allocatable :: Bex,Bey,Bez

  real(8), dimension(:), allocatable :: x,y,z

  integer :: i,j,k,nx,ny,nz
  real(8) :: hx,hy,hz

  ! Resistive intermediates
  real(8) :: Gax,Gay,Gaz,omega1,omega2,BLG
  real(8) :: cBfx,cBfy,cBfz
  real(8) :: rate_x_res,rate_y_res,rate_z_res

  ! Hall intermediates
  real(8) :: cJxB_x,cJxB_y,cJxB_z
  real(8) :: JxB_dot_b
  real(8) :: rate_x_hall,rate_y_hall,rate_z_hall

  ! File handling
  integer :: fileind, start_ind, end_ind
  character(4) :: str_ind, start_str, end_str

  call GET_COMMAND_ARGUMENT(1,start_str)
  call GET_COMMAND_ARGUMENT(2,end_str)
  read(start_str,'(I5)') start_ind
  read(end_str,'(I5)') end_ind

  do fileind = start_ind,end_ind

    write(str_ind,'(I4)') fileind
    print *, 'reading: output_'//trim(adjustl(str_ind))//'.nc'

    call load_field_from_netcdf('output_'//trim(adjustl(str_ind))//'.nc',x,y,z,Bx,By,Bz)

    nx=size(x); ny=size(y); nz=size(z)

    allocate(magB(nx,ny,nz), Fx(nx,ny,nz), Fy(nx,ny,nz), Fz(nx,ny,nz), magF(nx,ny,nz))
    allocate(Jx(nx,ny,nz), Jy(nx,ny,nz), Jz(nx,ny,nz))
    allocate(Bfx(nx,ny,nz), Bfy(nx,ny,nz), Bfz(nx,ny,nz))
    allocate(alpha(nx,ny,nz), lambda(nx,ny,nz))
    allocate(rate(nx,ny,nz))

    allocate(JxB_x(nx,ny,nz), JxB_y(nx,ny,nz), JxB_z(nx,ny,nz))
    allocate(Bex(nx,ny,nz), Bey(nx,ny,nz), Bez(nx,ny,nz))

    hx = abs(x(2)-x(1))
    hy = abs(y(2)-y(1))
    hz = abs(z(2)-z(1))

    !==============================
    ! First loop: base quantities
    !==============================
    do k=2,nz-1
    do j=2,ny-1
    do i=2,nx-1

      magB(i,j,k) = sqrt(Bx(i,j,k)**2 + By(i,j,k)**2 + Bz(i,j,k)**2)

      Jx(i,j,k) = 1/(2*hy)*(Bz(i,j+1,k)-Bz(i,j-1,k)) - 1/(2*hz)*(By(i,j,k+1)-By(i,j,k-1))
      Jy(i,j,k) = 1/(2*hz)*(Bx(i,j,k+1)-Bx(i,j,k-1)) - 1/(2*hx)*(Bz(i+1,j,k)-Bz(i-1,j,k))
      Jz(i,j,k) = 1/(2*hx)*(By(i+1,j,k)-By(i-1,j,k)) - 1/(2*hy)*(Bx(i,j+1,k)-Bx(i,j-1,k))

      Fx(i,j,k) = Jy(i,j,k)*Bz(i,j,k) - Jz(i,j,k)*By(i,j,k)
      Fy(i,j,k) = Jz(i,j,k)*Bx(i,j,k) - Jx(i,j,k)*Bz(i,j,k)
      Fz(i,j,k) = Jx(i,j,k)*By(i,j,k) - Jy(i,j,k)*Bx(i,j,k)
      magF(i,j,k) = sqrt(Fx(i,j,k)**2 + Fy(i,j,k)**2 + Fz(i,j,k)**2)

      Bfx(i,j,k) = (By(i,j,k)*Fz(i,j,k)-Fy(i,j,k)*Bz(i,j,k))/magF(i,j,k)
      Bfy(i,j,k) = (Bz(i,j,k)*Fx(i,j,k)-Fz(i,j,k)*Bx(i,j,k))/magF(i,j,k)
      Bfz(i,j,k) = (Bx(i,j,k)*Fy(i,j,k)-Fx(i,j,k)*By(i,j,k))/magF(i,j,k)

      alpha(i,j,k)  = (Jx(i,j,k)*Bx(i,j,k)+Jy(i,j,k)*By(i,j,k)+Jz(i,j,k)*Bz(i,j,k))/magB(i,j,k)**2
      lambda(i,j,k) = (Jx(i,j,k)*Bfx(i,j,k)+Jy(i,j,k)*Bfy(i,j,k)+Jz(i,j,k)*Bfz(i,j,k))/magB(i,j,k)**2

      JxB_x(i,j,k) = Fx(i,j,k)
      JxB_y(i,j,k) = Fy(i,j,k)
      JxB_z(i,j,k) = Fz(i,j,k)

      Bex(i,j,k) = Bx(i,j,k)/magB(i,j,k)
      Bey(i,j,k) = By(i,j,k)/magB(i,j,k)
      Bez(i,j,k) = Bz(i,j,k)/magB(i,j,k)

    end do
    end do
    end do

    !==============================
    ! Second loop: rate
    !==============================
    do k=3,nz-2
    do j=3,ny-2
    do i=3,nx-2

      ! --- Resistive ---
      cBfx = 1/(2*hy)*(Bfz(i,j+1,k)-Bfz(i,j-1,k)) - 1/(2*hz)*(Bfy(i,j,k+1)-Bfy(i,j,k-1))
      cBfy = 1/(2*hz)*(Bfx(i,j,k+1)-Bfx(i,j,k-1)) - 1/(2*hx)*(Bfz(i+1,j,k)-Bfz(i-1,j,k))
      cBfz = 1/(2*hx)*(Bfy(i+1,j,k)-Bfy(i-1,j,k)) - 1/(2*hy)*(Bfx(i,j+1,k)-Bfx(i,j-1,k))

      omega1 = (cBfx*Fx(i,j,k)+cBfy*Fy(i,j,k)+cBfz*Fz(i,j,k))/magF(i,j,k)
      omega2 = (cBfx*Bfx(i,j,k)+cBfy*Bfy(i,j,k)+cBfz*Bfz(i,j,k))/magB(i,j,k)**2

      BLG = (1/(2*hx)*(lambda(i+1,j,k)-lambda(i-1,j,k))*Bx(i,j,k) &
           + 1/(2*hy)*(lambda(i,j+1,k)-lambda(i,j-1,k))*By(i,j,k) &
           + 1/(2*hz)*(lambda(i,j,k+1)-lambda(i,j,k-1))*Bz(i,j,k))

      Gax = 1/(2*hx)*(alpha(i+1,j,k)-alpha(i-1,j,k))
      Gay = 1/(2*hy)*(alpha(i,j+1,k)-alpha(i,j-1,k))
      Gaz = 1/(2*hz)*(alpha(i,j,k+1)-alpha(i,j,k-1))

      rate_x_res = -1/magB(i,j,k)*( (lambda(i,j,k)*omega1-BLG)*Fx(i,j,k)/magF(i,j,k) &
                       + lambda(i,j,k)*(alpha(i,j,k)+omega2)*Bfx(i,j,k) &
                       + (Gay*Bz(i,j,k)-Gaz*By(i,j,k)) )

      rate_y_res = -1/magB(i,j,k)*( (lambda(i,j,k)*omega1-BLG)*Fy(i,j,k)/magF(i,j,k) &
                       + lambda(i,j,k)*(alpha(i,j,k)+omega2)*Bfy(i,j,k) &
                       + (Gaz*Bx(i,j,k)-Gax*Bz(i,j,k)) )

      rate_z_res = -1/magB(i,j,k)*( (lambda(i,j,k)*omega1-BLG)*Fz(i,j,k)/magF(i,j,k) &
                       + lambda(i,j,k)*(alpha(i,j,k)+omega2)*Bfz(i,j,k) &
                       + (Gax*By(i,j,k)-Gay*Bx(i,j,k)) )

      ! --- Hall ---
      cJxB_x = 1/(2*hy)*(JxB_z(i,j+1,k)-JxB_z(i,j-1,k)) - 1/(2*hz)*(JxB_y(i,j,k+1)-JxB_y(i,j,k-1))
      cJxB_y = 1/(2*hz)*(JxB_x(i,j,k+1)-JxB_x(i,j,k-1)) - 1/(2*hx)*(JxB_z(i+1,j,k)-JxB_z(i-1,j,k))
      cJxB_z = 1/(2*hx)*(JxB_y(i+1,j,k)-JxB_y(i-1,j,k)) - 1/(2*hy)*(JxB_x(i,j+1,k)-JxB_x(i,j-1,k))

      JxB_dot_b = Bex(i,j,k)*cJxB_x + Bey(i,j,k)*cJxB_y + Bez(i,j,k)*cJxB_z

      rate_x_hall = -(cJxB_x - JxB_dot_b*Bex(i,j,k))/magB(i,j,k)
      rate_y_hall = -(cJxB_y - JxB_dot_b*Bey(i,j,k))/magB(i,j,k)
      rate_z_hall = -(cJxB_z - JxB_dot_b*Bez(i,j,k))/magB(i,j,k)

      rate(i,j,k) = sqrt( (rate_x_res+rate_x_hall)**2 &
                         +(rate_y_res+rate_y_hall)**2 &
                         +(rate_z_res+rate_z_hall)**2 )

    end do
    end do
    end do

    print *, 'Writing rate_'//trim(adjustl(str_ind))//'.nc'
    call save_rate_to_netcdf(x,y,z,rate,'rate_'//trim(adjustl(str_ind))//'.nc')

    deallocate(x,y,z,Bx,By,Bz,magB,Fx,Fy,Fz,magF,Jx,Jy,Jz,Bfx,Bfy,Bfz,alpha,lambda,rate)
    deallocate(JxB_x,JxB_y,JxB_z,Bex,Bey,Bez)

  end do

end program main