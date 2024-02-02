! RRI_RivSlo
! river and slope interaction

subroutine funcrs( hr, hs )
use globals
implicit none

real(8) hr(ny, nx), hs(ny, nx)

real(8) hrs ! discharge amount from slope to river [m/s]
real(8) hr_top, hs_top, h1, h2
real(8) mu1, mu2, mu3
real(8) hr_new, len, b, ar ! add v1.4
integer i, j, k, count ! add v1.4 (k, count)

mu1 = (2.d0 / 3.d0) ** (3.d0 / 2.d0)
mu2 = 0.35d0
mu3 = 0.91d0

do i = 1, ny
 do j = 1, nx

  if( domain(i,j) .eq. 0 ) cycle
  if( riv(i,j) .eq. 0 ) cycle

  hs_top = hs(i,j)
  hr_top = hr(i,j) - depth(i,j)

  k = riv_ij2idx(i, j)
  len = len_riv_idx(k)

  ! (Case a) : (height = 0 and hr_top < 0) or (height > 0 and hr_top < 0 and hs_top <= height)
  ! -> From slope to river : step fall (hrs : positive)

  if( ( height(i,j) .eq. 0.d0 .and. hr_top .lt. 0.d0 ) .or.   &
      ( height(i,j) .gt. 0.d0 .and. hr_top .lt. 0.d0 .and. hs_top .le. height(i,j) ) ) then

   !hrs = mu1 * hs_top * sqrt( 9.81d0 * hs_top ) * dt * len / area
   hrs = mu1 * hs_top * sqrt( 9.81d0 * hs_top ) * dt * len * 2.d0 / area ! modified 1.4.2.4 
   if( hrs .gt. hs(i,j) ) hrs = hs(i,j)
   hs(i,j) = hs(i,j) - hrs
   !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
   call hr_update(hr(i, j), hrs * area, k, hr_new) ! add v1.4
   hr(i, j) = hr_new ! add v1.4
   qrs(i, j) = hrs

   ! avoid the situation of hr_top > hs_top
   hs_top = hs(i,j)
   hr_top = hr(i,j) - depth(i,j)
   if( hr_top .ge. -0.00001d0 .and. hr_top .gt. hs_top ) then
    do count = 1, 10
     call sec_h2b(hr(i,j), k, b)
     ar = len * b / area
     hrs = ( hs_top - hr_top ) / ( 1.d0 + 1.d0 / ar )
     hs(i,j) = hs(i,j) - hrs
     !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
     call hr_update(hr(i, j), hrs * area, k, hr_new) ! add v1.4
     hr(i, j) = hr_new ! add v1.4
     qrs(i, j) = qrs(i, j) + hrs
     if( abs(hs(i, j) - (hr(i,j) - depth(i, j))) .lt. 0.00001d0 ) exit
     hs_top = hs(i, j)
     hr_top = hr(i, j) - depth(i, j)
    enddo
    hr(i, j) = hs(i, j) + depth(i, j) ! 最終手段
   endif

  ! (Case b) : height > 0 and hr_top <= height and hr_top >= 0
  ! -> No exchange

  elseif( height(i,j) .gt. 0.d0 .and. hs_top .le. height(i,j) .and. hr_top .le. height(i,j) .and. hr_top .ge. 0.d0 ) then
   qrs(i, j) = 0.d0
  continue


  ! (Case c) : hs <= hrt & hrt >= height
  ! -> From river to slope : overtopping (hrs : negative)
  ! (incl. hs = 0 and hrt > 0)

  elseif( hs_top .le. hr_top .and. hr_top .ge. height(i,j) ) then

   h1 = hr_top - height(i, j)
   h2 = hs_top - height(i, j)
   if( h2 / h1 .le. 2.d0 / 3.d0 ) then
    !hrs = - mu2 * h1 * sqrt( 2.d0 * 9.81d0 * h1 ) * dt * len / area ! modified v1.4
    hrs = - mu2 * h1 * sqrt( 2.d0 * 9.81d0 * h1 ) * dt * len * 2.d0 / area ! modified v1.4.2.4
    else
    !hrs = - mu3 * h2 * sqrt( 2.d0 * 9.81d0 * (h1 - h2) ) * dt * len / area ! modified v1.4
    hrs = - mu3 * h2 * sqrt( 2.d0 * 9.81d0 * (h1 - h2) ) * dt * len * 2.d0 / area ! modified v1.4.2.4
   endif

   call sec_h2b(hr(i,j), k, b)
   ar = len * b / area
   if( abs(hrs / ar) .gt. (hr_top-height(i,j)) ) hrs = - (hr_top-height(i,j)) * ar
   !if( abs(hrs / area_ratio(i, j)) .gt. (hr_top-height(i,j)) ) hrs = - (hr_top-height(i,j)) * area_ratio(i, j)
   qrs(i, j) = hrs

   hs(i,j) = hs(i,j) - hrs
   !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
   call hr_update(hr(i, j), hrs * area, k, hr_new) ! add v1.4
   hr(i, j) = hr_new ! add v1.4

   ! avoid the situation of hs_top > hr_top
   hs_top = hs(i,j)
   hr_top = hr(i,j) - depth(i,j)
   !if( hr_top .gt. -0.00001d0 .and. hs_top .gt. hr_top ) then ! modified from v1.4 (ここはどっちか？)
   if( hs_top .gt. hr_top ) then
    do count = 1, 10
     call sec_h2b(hr(i,j), k, b)
     ar = len * b / area
     hrs = ( hs_top - hr_top ) / ( 1.d0 + 1.d0 / ar )
     hs(i,j) = hs(i,j) - hrs
     !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
     call hr_update(hr(i, j), hrs * area, k, hr_new) ! add v1.4
     hr(i, j) = hr_new ! add v1.4
     qrs(i, j) = qrs(i, j) + hrs
     if( abs(hs(i, j) - (hr(i,j) - depth(i, j))) .lt. 0.00001d0 ) exit
     hs_top = hs(i, j)
     hr_top = hr(i, j) - depth(i, j)
    enddo
    hr(i, j) = hs(i, j) + depth(i, j) ! 最終手段
   endif


  ! (Case d) : hs > hrt & hs >= height
  ! -> From slope to river : overtopping (hrs : positive)
  ! (incl. hs = 0 and hrt > 0)

  elseif( hs_top .ge. hr_top .and. hs_top .ge. height(i,j) ) then

   h1 = hs_top - height(i, j)
   h2 = hr_top - height(i, j)
   if( h2 / h1 .le. 2.d0 / 3.d0 ) then
    !hrs = mu2 * h1 * sqrt( 2.d0 * 9.81d0 * h1 ) * dt * len / area ! v1.4
    hrs = mu2 * h1 * sqrt( 2.d0 * 9.81d0 * h1 ) * dt * len * 2.d0 / area ! v1.4.2.4
   else
    !hrs = mu3 * h2 * sqrt( 2.d0 * 9.81d0 * (h1 - h2) ) * dt * len / area ! v1.4
    hrs = mu3 * h2 * sqrt( 2.d0 * 9.81d0 * (h1 - h2) ) * dt * len * 2.d0 / area ! v1.4.2.4
   endif

   if( hrs .gt. (hs_top - height(i,j)) ) hrs = hs(i,j) - height(i,j)
   qrs(i, j) = hrs

   hs(i,j) = hs(i,j) - hrs
   !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
   call hr_update(hr(i, j), hrs * area, k, hr_new) ! add v1.4
   hr(i, j) = hr_new ! add v1.4
   
   ! avoid the situation of hr_top > hs_top
   hs_top = hs(i,j)
   hr_top = hr(i,j) - depth(i,j)
   if( hr_top .ge. -0.00001d0 .and. hr_top .gt. hs_top ) then
    do count = 1, 10
     call sec_h2b(hr(i,j), k, b)
     ar = len * b / area
     hrs = ( hs_top - hr_top ) / ( 1.d0 + 1.d0 / ar )
     hs(i,j) = hs(i,j) - hrs
     !hr(i,j) = hr(i,j) + hrs / area_ratio(i,j)
     call hr_update(hr(i, j), hrs * area, k, hr_new) ! add v1.4
     hr(i, j) = hr_new ! add v1.4
     qrs(i, j) = qrs(i, j) + hrs
     if( abs(hs(i, j) - (hr(i,j) - depth(i, j))) .lt. 0.00001d0 ) exit
     hs_top = hs(i, j)
     hr_top = hr(i, j) - depth(i, j)
    enddo
    hr(i, j) = hs(i, j) + depth(i, j) ! 最終手段
   endif

  else

   ! Condition not considered above
   stop "Error : RivSlo"

  endif

  qrs(i, j) = qrs(i, j) / dble(dt) ! [m/s]

 enddo
enddo

end subroutine funcrs
