subroutine RRI_Div(qr_idx, hr_idx, qr_div_idx)
use globals
implicit none

real(8) qr_idx(riv_count), hr_idx(riv_count), qr_div_idx(riv_count)

integer l, k, kk, kk_div
real(8) zb_p, hr_p, zb_n, hr_n
real(8) dh, distance, qr_temp, hw
integer dif_p, dif_n


do l = 1, div_id_max
 k = div_org_idx(l)
 kk = down_riv_idx(k)
 kk_div = div_dest_idx(l)

 if( div_rate(l) .gt. 0.d0 ) then
  if( qr_idx(k) .gt. 0.d0 ) then
   qr_div_idx(k) = qr_idx(k) * div_rate(l)
   !qr_idx(k) = qr_idx(k) - qr_div_idx(k) ! comment out by T.Sayama on Dec 7, 2022 v.1.4.2.7
   !qr_idx(kk_div) = qr_idx(kk_div) + qr_div_idx(k)
  endif

 else ! div_rate(l) = 0.d0

  zb_p = zb_riv_idx(k)
  hr_p = hr_idx(k)
  dif_p = dif_riv_idx(k)

  distance = dis_riv_idx(k)

  ! information of the destination cell
  zb_n = zb_riv_idx(kk_div)
  hr_n = hr_idx(kk_div)
  dif_n = dif_riv_idx(kk_div)

  ! diffusion wave
  dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance ! diffussion

  ! kinematic wave
  if( dif_p .eq. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

  if( dh .ge. 0.d0 ) then
   hw = hr_p
   if( zb_p .lt. zb_n ) hw = max(0.d0, zb_p + hr_p - zb_n)
   call hq_riv(hw, dh, k, width_idx(k), qr_temp)
   qr_div_idx(k) = qr_temp
  else
   ! reverse flow
   hw = hr_n
   if( zb_n .lt. zb_p ) hw = max(0.d0, zb_n + hr_n - zb_p)
   dh = abs(dh)
   call hq_riv(hw, dh, kk_div, width_idx(k), qr_temp)
   qr_div_idx(k) = -qr_temp
  endif
  !qr_idx(k) = qr_idx(k) - qr_div_idx(k) ! comment out by T.Sayama on Dec 7, 2022 v.1.4.2.7
  !qr_idx(kk_div) = qr_idx(kk_div) + qr_div_idx(k)
 endif
enddo

end subroutine RRI_Div
