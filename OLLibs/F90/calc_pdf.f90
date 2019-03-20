!###############################################################
!
! The following routine takes an array as input and calculates the
! probability distribution functions with uniform statistical weight.
! This approach is diferent from computing histograms: uniform bin 
! size as the number of points in each bin might differ. This 
! routine makes the statistics in each bin equal.
!
!###############################################################
!
subroutine calc_pdf(a,nw,an,apdf,nx,ny,nz)
   Integer, Intent(In) :: nw,nx,ny,nz
   double precision, dimension(nx,ny,nz), intent(in) :: a 
   double precision, dimension(nx*ny*nz/nw), intent(out) :: an, apdf
!f2py    intent (in) nw,nb,nx,ny,nz
!f2py    intent (in) a
!f2py    intent (out) an,apdf
   Integer                           :: i,NNN,nb
   double precision                  :: rmsval
   double precision, dimension(nx*ny*nz) :: tmp
   double precision, dimension(nx,ny,nz) :: tmparr
   
   NNN=nx*ny*nz; nb=NNN/nw 
   apdf(:)=0.; an(:)=0.  ! Initialize the value array and pdf to 0.
   ! Take out the DC component and divide by the rms value of fluctuations
   tmparr=a-sum(a)/NNN; rmsval=dsqrt(sum(tmparr**2)/(nx*ny*nz))
   if (rmsval .ne. 0) tmparr=tmparr/rmsval
   ! Map the 2D array to 1D and sort it
   tmp=reshape(tmparr,(/NNN/))
   call qsort(tmp,NNN)
   ! Calculate the pdfs
   do i=1,nb
    an(i) = sum(tmp((i-1)*nw:(i-1)*nw+nw))/nw
    apdf(i) =nw/(maxval(tmp((i-1)*nw:(i-1)*nw+nw))- minval(tmp((i-1)*nw:(i-1)*nw+nw)))
   enddo
   apdf=apdf/NNN
End Subroutine calc_pdf
