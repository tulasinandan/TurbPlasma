subroutine cutit(ar,v,dt,lx,ly,lz,xi,yi,zi,theta,phi,nx,ny,nz,x,y,z,ao)
   implicit none
   integer, intent(in) :: nx,ny,nz
   double precision, intent(in) :: v, dt, lx, ly, lz, theta, phi,xi,yi,zi
   double precision, intent(in), dimension(nx,ny,nz) :: ar
  Real*8, Dimension(0:nx+1,0:ny+1,0:nz+1)       :: tmp
   double precision, parameter :: PI=4.D0*DATAN(1.D0)
   double precision :: dxg, dyg, dzg, xo, yo, zo, xn, yn, zn, dx, dy, &
                       dz, delta_l,st,ct,sp,cp
   double precision, intent(out), dimension(nx*ny*nz) :: x,y,z,ao
   integer :: i,nnn

! Copy guard-cells to tmp
   tmp(1:nx,1:ny,1:nz) = ar
   tmp(0,:,:)          = ar(nx,:,:)
   tmp(nx+1,:,:)       = ar(1,:,:)
   tmp(:,0,:)          = ar(:,ny,:)
   tmp(:,ny+1,:)       = ar(:,1,:)
   tmp(:,:,0)          = ar(:,:,nz)

   tmp(0,0,0)          = ar(nx,ny,nz)
   tmp(nx+1,0,0)       = ar(1,ny,nz)
   tmp(0,ny+1,0)       = ar(nx,1,nz)
   tmp(0,0,nz+1)       = ar(nx,ny,1)
   tmp(nx+1,ny+1,0)    = ar(1,1,nz)
   tmp(nx+1,0,nz+1)    = ar(1,ny,1)
   tmp(0,ny+1,nz+1)    = ar(nx,1,1)
   tmp(nx+1,ny+1,nz+1) = ar(1,1,1)

   delta_l= v*dt; nnn=nx*ny*nz
   st=sin(pi*theta/180.); ct=cos(pi*theta/180.)
   sp=sin(pi*phi/180.); cp=cos(pi*phi/180.)
   dx = delta_l*st*cp
   dy = delta_l*st*sp
   dz = delta_l*ct

  !write(*,*) dx, dy, dz

   dxg=lx/nx; dyg=ly/ny; dzg=lz/nz

   ao = 0.; i=1
   ! Initial x,y,z and the value of a_output at that point
   xo=xi; yo=yi; zo=zi
   x(i)=xi; y(i)=yi; z(i)=zi 
   call interp((/ xo, yo, zo /),nx, ny, nz,dxg,dyg,dzg,tmp,ao(i))

   do
     !if ((i .gt. 1.5*nx) .and. ((xn .lt. 0.1*lx) &
     !    .and. (yn .lt. 0.1*ly))) then
      if ((i .gt. nnn) .and. ((xn .gt. xi+0.1) &
          .and. (yn .gt. yi+0.1))) then
      !write(*,*) i,xn,yn
       exit
      elseif (i .gt. nnn) then
       write(*,*) 'Reached maximum number of points'
       exit
      endif

      xn = xo + dx; yn = yo + dy; zn = zo + dz 
      if (xn > lx) xn = xn - lx
      if (yn > ly) yn = yn - ly
      if (zn > lz) zn = zn - lz
!     write(*,*) xn,yn,zn
      call interp((/ xn, yn, zn /),nx, ny, nz,dxg,dyg,dzg,tmp,ao(i))
      xo=xn; yo=yn; zo=zn; i = i + 1
      x(i)=xo; y(i)=yo; z(i)=zo
   enddo

end subroutine cutit

!------------------------------------------------------------------------------
!             Interpolate quantities to particle positions
!------------------------------------------------------------------------------

! Given a position (x,y,z) and a quantity defined on grid points
! (e.g., jx), find the value of the quantity at that position via
! interpolation.
Subroutine interp(position,lnx,lny,lnz,dx,dy,dz,quantity,output)
  Implicit None
  Real*8, Dimension(3) :: position
  integer, intent(in) :: lnx, lny, lnz
  Real*8, intent(in) :: dx, dy, dz
  Real*8, Dimension(0:lnx+1,0:lny+1,0:lnz+1), intent(in) :: quantity
  Integer :: xval,yval,zval
  Real*8 :: xm,ym,zm,xp,yp,zp
  Real*8, intent(out) :: output

! Interpolate density and current.  First find gridpoint _below_ point.
  xval = int(position(1)/dx + 0.5)
  yval = int(position(2)/dy + 0.5)
  zval = int(position(3)/dz + 0.5)
! write(*,*) xval,yval,zval
  xm=position(1)/dx-xval+0.5; xp=1.-xm
  ym=position(2)/dy-yval+0.5; yp=1.-ym
  zm=position(3)/dz-zval+0.5; zp=1.-zm
! write(*,*) xm,xp, ym,yp,zm,zp

  output = (quantity(xval,yval,zval)          *xp*yp*zp +&
            quantity(xval,yval,zval+1)        *xp*yp*zm +&
            quantity(xval+1,yval,zval)        *xm*yp*zp +&
            quantity(xval+1,yval,zval+1)      *xm*yp*zm +&
            quantity(xval,yval+1,zval)        *xp*ym*zp +&
            quantity(xval,yval+1,zval+1)      *xp*ym*zm +&
            quantity(xval+1,yval+1,zval)      *xm*ym*zp +&
            quantity(xval+1,yval+1,zval+1)    *xm*ym*zm)

  return
End Subroutine interp

