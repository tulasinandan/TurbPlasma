! -*- f90 -*-
!#########################################
! SUBROUTINES TO COMPUTE CORRELATIONS 
!#########################################

subroutine correlation(a,b,axis,dx,step,rlen,nx,ny,nz,r,corr)
   implicit none
   integer, intent (in) :: nx,ny,nz,axis,rlen,step
   double precision, intent (in) :: dx
   double precision, intent (in), dimension(nx,ny,nz) :: a,b
   double precision, intent (out),dimension(rlen)     :: r,corr
!f2py intent (in) nx,ny,nz,axis,rlen,step
!f2py intent (in) dx
!f2py intent (in) a,b
!f2py intent (out) r,corr
   integer :: ii,x,y,z,xx,yy,zz
   double precision :: sumab
!  double precision, dimension(nx,ny,nz) :: br
   
   ! We can work with total sums instead of means as nx*ny*nz 
   ! cancels out in the final definition
!  sumab=sum(a*b)
    ! SELECT CASE (axis+1)
    !    CASE (1)
    !       print *, 'correlation along x'
    !    CASE (2)
    !       print *, 'correlation along y'
    !    CASE(3)
    !       print *, 'correlation along z'
    ! END SELECT

   do ii=0,rlen-1
   corr(ii+1) = 0.
   do z=1,nz; do y=1,ny; do x=1,nx
      SELECT CASE (axis+1)
         CASE (1)
            call roll_idx(x,nx,ii,xx); yy=y; zz=z
         CASE (2)
            xx=x; call roll_idx(y,ny,ii,yy); zz=z
         CASE(3)
            xx=x; yy=y; call roll_idx(z,nz,ii,zz)
      END SELECT
      corr(ii+1) = corr(ii+1)+a(x,y,z)*b(xx,yy,zz)
   enddo; enddo; enddo
!     call rollarr(b,axis+1,ii*step,br,nx,ny,nz)
!     corr(ii+1) = corr(ii+1)/sumab
      r(ii+1)    = ii*dx*step
   enddo
   corr(:)=corr(:)/corr(1)
end subroutine correlation
