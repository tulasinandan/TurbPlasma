!Routine to fill data gaps in a series
!Input:
!       time: the time series that is not linear in time
!       dt: expected time lag between succesive observations
!       tol: tolerance for dt
!       ai: A series to fill in the gaps alongside the time series
!       ln: length of ai
!       nmd: an approximation to extend the current length of time and ai by to
!       account for filling in the time gaps
!       fillval: the value to fill into ai when gaps are filled in time
!Output:
!       tn: filled time series
!       dn: filled ai series

subroutine fill_gaps(time, dt, tol, ai, ln, nmd, fillval, tn, dn,counter)
   double precision, intent(in) :: fillval, dt, tol
   integer, intent(in) :: ln, nmd
   integer*8, intent(out) :: counter
   double precision, intent(in), dimension(ln) :: ai, time
   double precision, intent(out), dimension(ln+nmd*2) :: tn, dn
   integer*8 :: i
   double precision :: del
   counter = 0
   do i=1,ln
      del = time(i+1) - time(i)
      if (abs(del-dt) <= tol) then
         dn(counter) = ai(i)
         tn(counter) = time(i)
         counter = counter + 1
      else
         do while (abs(time(i+1) - tn(counter-1)) > dt+tol)
            if (tn(counter-1)+dt >= time(i+1)) then
               EXIT
            else
               tn(counter) = tn(counter-1)+dt
               dn(counter) = fillval
               counter = counter + 1
            endif
         enddo
      endif
   enddo
end subroutine fill_gaps
