! -*- f90 -*-
subroutine strfn_vec(ax,ay,az,nt,lags,nlag,orders,nord,dt,r,S)
   implicit none
   integer*8, intent(in) :: nt,nlag,nord
   double precision, intent(in), dimension(nt) :: ax,ay,az
   integer*8, intent(in), dimension(nlag) :: lags
   integer*8, intent(in), dimension(nord) :: orders
   double precision, intent(in) :: dt
   integer*8 :: i,k,l,n
   double precision, intent(out), dimension(nord,nlag) :: S
   double precision, intent(out), dimension(nlag) :: r
   double precision :: dux, duy, duz, cnt, absu

   S = 0.0
   do k=1,nlag
      l=lags(k)
      cnt=nt-l
      do i=l+1,nt
         if ( isnan(ax(i)) .or. isnan(ax(i-l)) &
         .or. isnan(ay(i)) .or. isnan(ay(i-l)) &
         .or. isnan(az(i)) .or. isnan(az(i-l)) &
            ) then
            cnt=cnt-1
         else
            dux = ax(i)-ax(i-l)
            duy = ay(i)-ay(i-l)
            duz = az(i)-az(i-l)
            absu = sqrt(dux**2+duy**2+duz**2)
         endif
         do n=1,nord
            S(n,k) = S(n,k)+absu**orders(n)
         enddo
      enddo
      S(:,k) = S(:,k)/cnt!(nt-l-1)
      r(k)   = l*dt
   enddo
end subroutine strfn_vec

subroutine strfn(ax,nt,lags,nlag,orders,nord,dt,r,S)
   implicit none
   integer*8, intent(in) :: nt,nlag,nord
   double precision, intent(in), dimension(nt) :: ax
   integer*8, intent(in), dimension(nlag) :: lags
   integer*8, intent(in), dimension(nord) :: orders
   double precision, intent(in) :: dt
   integer*8 :: i,k,l,n
   double precision, intent(out), dimension(nord,nlag) :: S
   double precision, intent(out), dimension(nlag) :: r
   double precision :: cnt, absu

   S = 0.0
   do k=1,nlag
      l=lags(k)
      cnt=nt-l
      do i=l+1,nt
         if (isnan(ax(i)) .or. isnan(ax(i-l))) then
            cnt=cnt-1
         else
            absu = abs(ax(i)-ax(i-l))
         endif
         do n=1,nord
            S(n,k) = S(n,k)+absu**orders(n)
         enddo
      enddo
      S(:,k) = S(:,k)/cnt!(nt-l-1)
      r(k)   = l*dt
   enddo
end subroutine strfn
! -*- f90 -*-
! Correlation func using Blackman-Tukey 

! Nested loop for finding the diagonal correlation functions
! x, y are the two (scalar) variables as function of time

subroutine corr_loop(nlag,M,x,y,R)  

 implicit none
 integer*8, intent(in) :: M,nlag 
 real(8) :: my_sum
 real(8), intent(in), dimension(0:M-1) :: x, y
 real(8), intent(out), dimension(0:nlag-1) :: R
 integer :: P, n 

 DO n = 0,nlag-1,1
   my_sum = 0.0
    DO P = 0,M-n-1,1
     my_sum = my_sum + x(P)*y(P+n)
    END DO
    R(n) = my_sum/real(M-n)
  END DO

end subroutine

subroutine correlation(a,b,nt,lags,nlag,dt,r,cr)
  !use ieee_arithmetic
   implicit none
   integer*8, intent(in) :: nt, nlag
   double precision, intent(in) :: dt
   integer*8, intent(in), dimension(nlag) :: lags
   double precision, intent(in), dimension(nt) :: a,b
   double precision, intent(out), dimension(nlag) :: r,cr
   double precision :: aav, bav
   integer*8 :: i,k,l,cnt
   
   do k=1,nlag
      l=lags(k)
      cr(k)=0; aav=0; bav=0; cnt=nt-l
      do i=1,nt-l
        !if (ieee_is_nan(a(i)) .and. ieee_is_nan(b(i+l)) then
         if (isnan(a(i)) .or. isnan(b(i+l))) then
            cnt=cnt-1
           !continue
         else
            cr(k) = cr(k)+a(i)*b(i+l)
            aav = aav+a(i)
            bav = bav+b(i+l)
!           cnt=cnt+1
         endif
      enddo
      r(k)  = l*dt
      cr(k) = (cr(k)*cnt - aav*bav)/cnt**2
   enddo
!  cr = cr/cr(1)
end subroutine correlation
