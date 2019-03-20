double precision Function nansum(a,n)
   implicit none
   integer*8, intent(in) :: n
   double precision, intent(in), dimension(n) :: a
   integer*8 :: i
   nansum=0.
   do i = 1,n
      if ( isnan(a(i)) ) cycle
      nansum=nansum+a(i)
   enddo
end Function nansum

double precision Function nanmean(a,n)
   implicit none
   integer*8, intent(in) :: n
   double precision, intent(in), dimension(n) :: a
   integer*8 :: i,counter
   nanmean=0.; counter=0
   do i = 1,n
      if ( isnan(a(i)) ) cycle
      nanmean=nanmean+a(i)
      counter=counter+1
   enddo
   if (counter >= 1) nanmean=nanmean/counter
end Function nanmean

double precision Function nanstd(a,n)
   implicit none
   integer*8, intent(in) :: n
   double precision, intent(in), dimension(n) :: a
   integer*8 :: i,counter
   double precision :: nmean,nanmean
   nmean=nanmean(a,n)
   nanstd=0.; counter=0
   do i=1,n
      if ( isnan(a(i)) ) cycle
      nanstd=nanstd+(a(i)-nmean)**2
      counter=counter+1
   enddo
   if (counter >= 1) nanstd=sqrt(nanstd/counter)
end Function nanstd

subroutine dropna(ai,n,ao,counter)
   implicit none
   integer*8, intent(in) :: n
   double precision, intent(in), dimension(n) :: ai
   double precision, intent(out), dimension(n) :: ao
   integer*8, intent(out)  :: counter
   integer*8 :: i
   ao=0.; counter=0
   do i=1,n
      if ( isnan(ai(i)) ) cycle
      counter = counter+1
      ao(counter) = ai(i)
   enddo
end subroutine dropna

double precision Function nan_mean(a,n)
   implicit none
   integer*8, intent(in) :: n
   double precision, intent(in), dimension(n) :: a
   double precision :: at(n)
   integer*8 :: cntr
   call dropna(a,n,at,cntr)
   nan_mean=sum(at(1:cntr))/cntr
end Function nan_mean

double precision Function nan_median(a,n)
   implicit none
   integer*8, intent(in) :: n
   double precision, intent(in), dimension(n) :: a
   double precision :: at(n)
   double precision, allocatable, dimension(:) :: tmp
   integer*8 :: cntr
   call dropna(a,n,at,cntr)
   allocate(tmp(cntr)); tmp=at(1:cntr)
   call qsort(tmp,cntr)
   nan_median=tmp(int(cntr/2+1))
   deallocate(tmp)
end Function nan_median
