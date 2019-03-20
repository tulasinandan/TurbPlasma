! -*- f90 -*-

!#############################################################
! COMPUTE 2D SPECTRUM OF AN ARRAY
!#############################################################

subroutine pspec2d(a,kp,kspec,axis,nk1,nk2,nx,ny,nz)
   integer, intent (in) :: nx,ny,nz,nk1,nk2,axis
   double precision, dimension(nx,ny,nz), intent (in) :: a
   double precision, dimension(nk1,nk2), intent (out) :: kp,kspec
!f2py intent (in) nx,ny,nz,axis,nk
!f2py intent (in) a
!f2py intent (out) kp, kspec
   ! LOCAL VARIABLES
   double complex  , allocatable, dimension(:,:,:) :: ffld
   double precision, allocatable, dimension(:,:,:) :: fenergy
   double precision, allocatable, dimension(:,:) :: kx, ky, kz
   integer*8 plan, LX1, NORM, i, j, k, ntop, nnx, nny

   include 'fftw3.f'

   LX1=nx/2 + 1; NORM = nx*ny*nz
   select case (axis+1)
      case (1)
        !allocate(fenp(ny,nz), kp(ny,nz), ky(ny,nz), kz(ny,nz))
         allocate(ky(ny,nz), kz(ny,nz))
         do i=1,ny; ky(i,:)=dble(mod(i-1+ny/2,ny)-ny/2); end do
         do i=1,nz; kz(:,i)=dble(mod(i-1+nz/2,nz)-nz/2); end do
         kp=sqrt(ky*ky+kz*kz); kp(1,1) = dble(1.e-6)
         ntop=min(ny,nz)/2; nnx=ny; nny=nz
         deallocate(ky,kz)
      case (2)
        !allocate(fenp(LX1,nz),kp(LX1,nz), kx(LX1,nz), kz(LX1,nz))
         allocate(kx(LX1,nz), kz(LX1,nz))
         do i=1,LX1; kx(i,:)=dble(i-1.0); enddo
         do i=1,nz;  kz(:,i)=dble(mod(i-1+nz/2,nz)-nz/2); end do
         kp=sqrt(kx*kx+kz*kz); kp(1,1) = dble(1.e-6)
         ntop=min(nx,nz)/2; nnx=LX1; nny=nz
         deallocate(kx,kz)
      case (3)
        !allocate(fenp(LX1,ny),kp(LX1,ny), kx(LX1,ny), ky(LX1,ny))
         allocate(kx(LX1,ny), ky(LX1,ny))
         do i=1,LX1; kx(i,:)=dble(i-1.0); enddo
         do i=1,ny;  ky(:,i)=dble(mod(i-1+ny/2,ny)-ny/2); end do
         kp=sqrt(kx*kx+ky*ky); kp(1,1) = dble(1.e-6)
         ntop=min(nx,ny)/2; nnx=LX1; nny=ny
         deallocate(kx,ky)
   end select


   allocate(ffld(LX1,ny,nz), fenergy(LX1,ny,nz))
   ! CREATE THE FFTW PLAN
   call dfftw_plan_dft_r2c_3d(plan,Nx,Ny,Nz,a,ffld,FFTW_ESTIMATE)
   call dfftw_execute_dft_r2c(plan, a, ffld)
   ffld=ffld/NORM ! normalize the quantities in fourier space

   ! Calculate the energy in the fourier space, make sure 
   ! that (1,:) part of the array is multiplied by half.
   fenergy=ffld*conjg(ffld); fenergy(1,:,:)=0.5*fenergy(1,:,:)
   kspec = sum(fenergy,dim=axis+1)
   ! DESTROY THE PLAN AND DEALLOCATE ARRAYS
   call dfftw_destroy_plan(plan)
   deallocate(ffld,fenergy)
end subroutine pspec2d

!#############################################################
! COMPUTE SPECTRUM OF AN ARRAY
!#############################################################

subroutine perpspec(a,kk,kspec,axis,nk,nx,ny,nz)
   integer, intent (in) :: nx,ny,nz,nk,axis
   double precision, dimension(nx,ny,nz), intent (in) :: a
   double precision, dimension(0:nk), intent (out) :: kk,kspec
!f2py intent (in) nx,ny,nz,axis,nk
!f2py intent (in) a
!f2py intent (out) kk, kspec
   ! LOCAL VARIABLES
   double complex  , allocatable, dimension(:,:,:) :: ffld
   double precision, allocatable, dimension(:,:,:) :: fenergy
   double precision, allocatable, dimension(:,:) :: fenp, kx, ky, kz, kp
   integer*8 plan, LX1, NORM, i, j, k, ntop, nnx, nny

   include 'fftw3.f'

   LX1=nx/2 + 1; NORM = nx*ny*nz
   select case (axis+1)
      case (1)
         allocate(fenp(ny,nz), kp(ny,nz), ky(ny,nz), kz(ny,nz))
         do i=1,ny; ky(i,:)=dble(mod(i-1+ny/2,ny)-ny/2); end do
         do i=1,nz; kz(:,i)=dble(mod(i-1+nz/2,nz)-nz/2); end do
         kp=sqrt(ky*ky+kz*kz); kp(1,1) = dble(1.e-6)
         ntop=min(ny,nz)/2; nnx=ny; nny=nz
         deallocate(ky,kz)
      case (2)
         allocate(fenp(LX1,nz),kp(LX1,nz), kx(LX1,nz), kz(LX1,nz))
         do i=1,LX1; kx(i,:)=dble(i-1.0); enddo
         do i=1,nz;  kz(:,i)=dble(mod(i-1+nz/2,nz)-nz/2); end do
         kp=sqrt(kx*kx+kz*kz); kp(1,1) = dble(1.e-6)
         ntop=min(nx,nz)/2; nnx=LX1; nny=nz
         deallocate(kx,kz)
      case (3)
         allocate(fenp(LX1,ny),kp(LX1,ny), kx(LX1,ny), ky(LX1,ny))
         do i=1,LX1; kx(i,:)=dble(i-1.0); enddo
         do i=1,ny;  ky(:,i)=dble(mod(i-1+ny/2,ny)-ny/2); end do
         kp=sqrt(kx*kx+ky*ky); kp(1,1) = dble(1.e-6)
         ntop=min(nx,ny)/2; nnx=LX1; nny=ny
         deallocate(kx,ky)
   end select


   allocate(ffld(LX1,ny,nz), fenergy(LX1,ny,nz))
   ! CREATE THE FFTW PLAN
   call dfftw_plan_dft_r2c_3d(plan,Nx,Ny,Nz,a,ffld,FFTW_ESTIMATE)
   call dfftw_execute_dft_r2c(plan, a, ffld)
   ffld=ffld/NORM ! normalize the quantities in fourier space

   ! Calculate the energy in the fourier space, make sure 
   ! that (1,:) part of the array is multiplied by half.
   fenergy=ffld*conjg(ffld); fenergy(1,:,:)=0.5*fenergy(1,:,:)
   fenp = sum(fenergy,dim=axis+1)
   kspec=0. ! Initialize the spectrum to zero
   ! Do the sum of energy over shells of k and save in the spectrum
   do i=1,nnx; do j=1,nny
      k=nint(kp(i,j))
      if (k <= ntop) then
         kspec(k)=kspec(k)+fenp(i,j)
         kk(k) = k
      endif
   enddo; enddo   
   ! DESTROY THE PLAN AND DEALLOCATE ARRAYS
   call dfftw_destroy_plan(plan)
   deallocate(ffld,fenergy,fenp,kp)
end subroutine perpspec

!###############################################################
! COMPUTE SPECTRUM OF AN ARRAY FOR A BOX WITH DIFFERENT LENGTHS
!###############################################################

subroutine lperpspec(a,axis,nk,lx,ly,lz,nx,ny,nz,kk,kspec)
   double precision, parameter :: pi=3.1415926535897932385
   integer, intent (in) :: nx,ny,nz,nk,axis
   double precision, intent (in) :: lx, ly, lz
   double precision, dimension(nx,ny,nz), intent (in) :: a
   double precision, dimension(0:nk), intent (out) :: kk,kspec
   ! LOCAL VARIABLES
   double complex  , allocatable, dimension(:,:,:) :: ffld
   double precision, allocatable, dimension(:,:,:) :: fenergy
   double precision, allocatable, dimension(:,:) :: fenp
   double precision, allocatable, dimension(:) :: kkx, kky
   double precision :: kminx, kminy, kminz, kp, kmin
   integer(kind=8) plan, FFTW_ESTIMATE, LX1, NORM, i, j, k, ktop, nnx, nny
   PARAMETER (FFTW_ESTIMATE=64)
!f2py intent (in) axis,nk,lx,ly,lz,nx,ny,nz
!f2py intent (in) a
!f2py intent (out) kk, kspec

   LX1=nx/2 + 1; NORM = nx*ny*nz
   kminx=2*pi/lx; kminy = 2*pi/ly; kminz = 2*pi/lz
   select case (axis+1)
      case (1)
         allocate(fenp(ny,nz), kkx(ny), kky(nz))
         !! kkx = ky, kky = kz
         do i=1,ny; kkx(i)=dble(mod(i-1+ny/2,ny)-ny/2)*kminy; end do
         do i=1,nz; kky(i)=dble(mod(i-1+nz/2,nz)-nz/2)*kminz; end do
         ktop=min(ny*kminy/2,nz*kminz/2); nnx=ny; nny=nz
         kmin=min(kminy,kminz)
      case (2)
         allocate(fenp(LX1,nz), kkx(LX1), kky(nz))
         !! kkx = kx, kky = kz
         do i=1,LX1; kkx(i)=dble(i-1.0)*kminx; enddo
         do i=1,nz;  kky(i)=dble(mod(i-1+nz/2,nz)-nz/2)*kminz; end do
         ktop=min(nx*kminx/2,nz*kminz/2); nnx=LX1; nny=nz
         kmin=min(kminx,kminz)
      case (3)
         allocate(fenp(LX1,ny), kkx(LX1), kky(ny))
         !! kkx = kx, kky = ky
         do i=1,LX1; kkx(i)=dble(i-1.0)*kminx; enddo
         do i=1,ny;  kky(i)=dble(mod(i-1+ny/2,ny)-ny/2)*kminy; end do
         ktop=min(nx*kminx/2,ny*kminy/2); nnx=LX1; nny=ny
         kmin=min(kminx,kminy)
   end select

   allocate(ffld(LX1,ny,nz), fenergy(LX1,ny,nz))
   ! CREATE THE FFTW PLAN
   call dfftw_plan_dft_r2c_3d(plan,nx,ny,nz,a,ffld,FFTW_ESTIMATE)
   call dfftw_execute_dft_r2c(plan, a, ffld)
   ffld=ffld/NORM ! normalize the quantities in fourier space

   ! Calculate the energy in the fourier space, make sure 
   ! that (1,:) part of the array is multiplied by half.
   fenergy=ffld*conjg(ffld); fenergy(1,:,:)=0.5*fenergy(1,:,:)
   fenp = sum(fenergy,dim=axis+1)
   kspec=0. ! Initialize the spectrum to zero
   ! Do the sum of energy over shells of k and save in the spectrum
   do i=1,nnx; do j=1,nny
      kp=dsqrt(kkx(i)**2+kky(j)**2)
      if (kp <= ktop) then
         k = int(dsqrt((kkx(i)/kmin)**2+(kky(j)/kmin)**2)+0.5)
         kspec(k)=kspec(k)+fenp(i,j)
         kk(k) = kp
      endif
   enddo; enddo   
   ! DESTROY THE PLAN AND DEALLOCATE ARRAYS
   call dfftw_destroy_plan(plan)
   deallocate(ffld,fenergy,fenp,kkx,kky)
end subroutine lperpspec
