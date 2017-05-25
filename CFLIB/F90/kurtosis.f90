! -*- f90 -*-

!###############################################################
!!   SUBROUTINE TO COMPUTE THE KURTOSIS OF AN ARRAY
!###############################################################
subroutine kurtosis(a,nx,ny,nz,krt)
   integer, intent(in) :: nx, ny, nz
   double precision, intent(out):: krt
   double precision, dimension(nx,ny,nz), intent(in) :: a
!f2py intent (in) nx,ny,nz
!f2py intent (in) a
!f2py intent (out) krt
   double precision, dimension(nx,ny,nz) :: temp
   double precision :: tmp, NORM
   NORM=nx*ny*nz
   tmp=sum(a)/NORM; temp=a-tmp
   tmp=sum(temp**2)/NORM
   krt=sum(temp**4)/(NORM*tmp**2)
end subroutine kurtosis

!###############################################################
!!   SUBROUTINE TO COMPUTE SCALE DEPENDENT KURTOSIS OF AN ARRAY
!###############################################################
subroutine sdk(a,axis,nsd,nx,ny,nz,sdkrt)
   integer, intent(in) :: nx,ny,nz,axis,nsd
   double precision, dimension(nx,ny,nz), intent(in) :: a
   double precision, dimension(nsd), intent(out) ::sdkrt
!f2py intent (in) a
!f2py intent (in) nx, ny, nz, axis, nsd
!f2py intent (out) sdkrt
   double precision, dimension(nx,ny,nz) :: tmp
   integer :: ii
   double precision :: rmsval
!
   do ii=1,nsd
      call rollarr(a,axis,ii,tmp,nx,ny,nz)
      rmsval=dsqrt(sum(tmp**2)/(nx*ny*nz))
      if (rmsval .ne. 0) tmp=tmp/rmsval
      call kurtosis(tmp,nx,ny,nz,sdkrt(ii))
   enddo
end subroutine sdk

!###############################################################
!!   FORTRAN EQUIVALENT OF CSHIFT BUT FASTER AS IT DOES NOT
!!   CREATE EXTRA COPIES
!###############################################################
subroutine rollarr(a,axis,inc,ar,nx,ny,nz)
   double precision, dimension(nx,ny,nz), intent(in) :: a
   integer, intent(in) :: nx,ny,nz,axis,inc
   double precision, dimension(nx,ny,nz), intent(out):: ar
!f2py intent (in) a
!f2py intent (in) nx,ny,nz,axis,inc
!f2py intent (out) ar
   integer :: x,y,z,xx,yy,zz
!
   do z=1,nz; do y=1,ny; do x=1,nx
      SELECT CASE (axis+1)
         CASE (1)
            call roll_idx(x,nx,inc,xx); yy=y; zz=z
         CASE (2)
            xx=x; call roll_idx(y,ny,inc,yy); zz=z
         CASE(3)
            xx=x; yy=y; call roll_idx(z,nz,inc,zz)
      END SELECT
      ar(x,y,z) = a(xx,yy,zz) 
   enddo; enddo; enddo
end subroutine rollarr

!###############################################################
!!   SUBROUTINE TO ROLL AN INDEX ON PERIODIC GRID
!###############################################################
subroutine roll_idx(x,nx,inc,xx)
   integer, intent(in) :: x,nx,inc
   integer, intent(out):: xx
!f2py intent (in) x,nx,inc
!f2py intent (out) xx
   xx = merge(x+inc, x-nx+inc+1, (x .lt. nx-inc))
!  if (x .lt. nx-inc) then 
!     xx=x+inc 
!  else 
!     xx=x-nx+inc+1
!  endif
end subroutine roll_idx

