subroutine vdist3d(vx,vy,vz,lenv,nx,ny,nz,dfn)
   implicit none
   integer, intent(in) :: lenv, nx, ny, nz
   double precision, intent(in), dimension(lenv) :: vx, vy, vz
   double precision, intent(out), dimension(nx,ny,nz) :: dfn
   integer :: i,j,k,l
   double precision :: dvx, dvy, dvz
   double precision, dimension(6) :: mm

   mm=(/ maxval(vx) , minval(vx) &
        ,maxval(vy) , minval(vy) &
        ,maxval(vz) , minval(vz) /)

   dvx = (mm(1) - mm(2))/nx
   dvy = (mm(3) - mm(4))/ny
   dvz = (mm(5) - mm(6))/nz

   write(*,*) dvx,dvy,dvz

   do i=1,lenv
      j=floor((vx(i)-mm(2))/dvx + 0.5)
      k=floor((vy(i)-mm(4))/dvy + 0.5)
      l=floor((vz(i)-mm(6))/dvz + 0.5)
      !write(*,*) 'j,k,l',j,k,l
      dfn(j,k,l) = dfn(j,k,l) + 1
   enddo
end subroutine vdist3d
