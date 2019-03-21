!Slide a window by one point to determine outliers/bad data
!Input:
!       ai: array in to remove bad data from
!       ws: Window size to slide through ai
!       sigma: the factor of sigma to remove outliers outside of the range
!       fill: the fill value to replace outlier values with
!             special value fill=-2.379152e33 fills with median value
!Output:
!       ao: array that has all outliers removed from ai

subroutine hampel(ai, ao, nt, ws, numdev, fill)
   implicit none
   real, intent(in) :: fill
   integer*8, intent(in) :: nt, ws, numdev
   double precision, intent(in), dimension(nt) :: ai
   double precision, intent(out), dimension(nt) :: ao
   double precision :: X0, tmp(2*ws+1), L, S0, nan_median
   integer*8 :: i, cntr, win
   L=1.4286
   win=2*ws+1
   ao=ai
   do i=1+ws,nt-ws
      call dropna(ai(i-ws:i+ws),win,tmp,cntr)
      if (cntr .eq. 0) then
         cycle
      else 
         X0 = nan_median(ai(i-ws:i+ws),win)
         S0 = L*nan_median(abs(ai(i-ws:i+ws)-X0),win)
         if (abs(ai(i)-X0) .gt. numdev*S0) then
            if (fill .eq. -2.379152e33) then
               ao(i) = X0
            else
               ao(i) = fill
            endif
         endif
      endif
   enddo
end subroutine
