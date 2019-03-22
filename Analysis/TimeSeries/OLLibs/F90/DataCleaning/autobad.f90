!Slide a window by one point to determine outliers/bad data
!Input:
!       ai: array in to remove bad data from
!       ws: Window size to slide through ai
!       sigma: the factor of sigma to remove outliers outside of the range
!       fill: the fill value to replace outlier values with
!             special value fill=-2.379152e33 fills with median value
!Output:
!       ao: array that has all outliers removed from ai

subroutine autobad(ai, ao, nt, ws, numdev, fill)
   implicit none
   real, intent(in) :: fill
   integer*8, intent(in) :: nt, ws, numdev
   double precision, intent(in), dimension(nt) :: ai
   double precision, intent(out), dimension(nt) :: ao
   double precision :: nanmean, nanstd, X0, S0
   integer*8 :: i, win
   win=2*ws+1
   ao=ai
   do i=1+ws,nt-ws
      X0 = nanmean(ao(i-ws:i+ws),win)
      S0 =  nanstd(ao(i-ws:i+ws),win)
      if (abs(ai(i)-X0) > numdev*S0) then
         if (fill .eq. -2.379152e33) then
            ao(i) = X0
         else
            ao(i) = fill
         endif
      endif
   enddo
end subroutine
