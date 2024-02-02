! RRI_Riv.f90

! variable definition (river)

subroutine funcr( vr_idx, fr_idx, qr_idx )
use globals
use dam_mod, only: dam_switch, damflg, dam_qout, dam_volmax, dam_floodq, dam_num, dam_loc, dam_qin
implicit none

real(8) hr_idx(riv_count), vr_idx(riv_count), fr_idx(riv_count), qr_idx(riv_count)

real(8) qr_sum_idx(riv_count), qr_div_idx(riv_count)
integer i, k, kk, itemp, jtemp , l
integer kk_div

fr_idx(:) = 0.d0
qr_idx(:) = 0.d0
qr_sum_idx(:) = 0.d0
qr_div_idx(:) = 0.d0

! add v1.4
do k = 1, riv_count
 call vr2hr(vr_idx(k), k, hr_idx(k)) 
enddo

! boundary condition for river (water depth boundary)
if( bound_riv_wlev_switch .ge. 1 ) then
 itemp = -1
 do jtemp = 1, tt_max_bound_riv_wlev
  if( t_bound_riv_wlev(jtemp-1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_riv_wlev(jtemp) ) itemp = jtemp
 enddo
 do k = 1, riv_count
  if( bound_riv_wlev_idx(itemp, k) .le. -100.0 ) cycle ! not boundary
  hr_idx(k) = bound_riv_wlev_idx(itemp, k)
  call hr2vr(hr_idx(k), k, vr_idx(k)) ! add v1.4
 enddo
endif

call qr_calc(hr_idx, qr_idx)
if( div_switch .eq. 1 ) call RRI_Div(qr_idx, hr_idx, qr_div_idx)

! dam control
if( dam_switch .eq. 1 ) then
 call dam_prepare(qr_idx) ! calculate inflow to dam
 do i = 1, dam_num
  if( dam_volmax(i) .gt. 0.d0 ) then
   ! dam
   call dam_operation( dam_loc(i) )
   qr_idx( dam_loc(i) ) = dam_qout(i)
   qr_sum_idx( dam_loc(i) ) = qr_sum_idx( dam_loc(i) ) + dam_qin( dam_loc(i) ) - dam_qout(i)
  elseif( dam_volmax(i) .eq. 0.d0 ) then
   ! barrage
   if( hr_idx( dam_loc(i) ) .le. dam_floodq(i) ) then
    qr_idx( dam_loc(i) ) = 0.d0
   endif
  else
   ! water gate (dam_volmax(i) < 0) ! added on Aug 7, 2021 by T.Sayama
   k = dam_loc(i)
   kk = down_riv_idx(k)
   call gate_operation( dam_loc(i), hr_idx(k), hr_idx(kk) )
   qr_idx( dam_loc(i) ) = dam_qout(i)
   !qr_sum_idx( dam_loc(i) ) = qr_sum_idx( dam_loc(i) ) + dam_qin( dam_loc(i) ) - dam_qout(i)
  endif
 enddo
endif

! boundary condition for river (discharge boundary)
if( bound_riv_disc_switch .ge. 1 ) then
 do jtemp = 1, tt_max_bound_riv_disc
  if( t_bound_riv_disc(jtemp-1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_riv_disc(jtemp) ) itemp = jtemp
 enddo
 do k = 1, riv_count
  if( bound_riv_disc_idx(itemp, k) .le. -100.0 ) cycle ! not boundary
  !qr_idx(k) = bound_riv_disc_idx(itemp, k) / area  ! qr_idx: discharge per unit area
  qr_idx(k) = bound_riv_disc_idx(itemp, k) ! modified v1.4
  ! linear interpolation of the boundary condition
  !qr_idx(k) = bound_riv_disc_idx(itemp-1, k) * (t_bound_riv_disc(itemp) - (time + ddt)) &
  !           + bound_riv_disc_idx(itemp, k) * ((time + ddt) - t_bound_riv_disc(itemp-1))
  !qr_idx(k) = qr_idx(k) / (t_bound_riv_disc(itemp) - t_bound_riv_disc(itemp-1))
  hr_idx(k) = 0.d0
  vr_idx(k) = 0.d0 ! add v1.4
 enddo
endif

! qr_sum > 0 --> discharge flowing out from a cell
do k = 1, riv_count
 ! outflow from (k)
 qr_sum_idx(k) = qr_sum_idx(k) + qr_idx(k)
 kk = down_riv_idx(k)
 if(domain_riv_idx(kk).eq.0) cycle
 ! qr_sum minus (flowing into) discharge at the destination cell
 qr_sum_idx(kk) = qr_sum_idx(kk) - qr_idx(k)
enddo

! diversion
if( div_switch .eq. 1 ) then
 do l = 1, div_id_max
 ! outflow from (k)
  k = div_org_idx(l)
  kk = down_riv_idx(k)
  if( div_rate(l) .le. 0.d0 ) then ! modified by T.Sayama on Dec. 7, 2022 v1.4.2.7
   qr_sum_idx(k) = qr_sum_idx(k) + qr_div_idx(k)
  else
   qr_sum_idx(kk) = qr_sum_idx(kk) + qr_div_idx(k)
  endif
  kk = div_dest_idx(l)
  if(domain_riv_idx(kk).eq.0) cycle
  ! qr_sum minus (flowing into) discharge at the destination cell
  qr_sum_idx(kk) = qr_sum_idx(kk) - qr_div_idx(k) 
 enddo
endif

! qr_sum divide by area_ratio
! (area_ratio = ratio of river area against total cell area) 
!fr_idx = -qr_sum_idx / area_ratio_idx
fr_idx = - qr_sum_idx ! modified v1.4

end subroutine funcr


! lateral discharge (river)

subroutine qr_calc(hr_idx, qr_idx)
use globals
use dam_mod, only: damflg, dam_floodq
implicit none

real(8) hr_idx(riv_count), qr_idx(riv_count), qr_div_idx(riv_count)

integer k, kk
real(8) zb_p, hr_p
real(8) zb_n, hr_n
real(8) dh, distance
real(8) qr_temp, hw
integer dif_p, dif_n

qr_idx(:) = 0.d0
qr_div_idx(:) = 0.d0

!$omp parallel do private(kk,zb_p,hr_p,distance,zb_n,hr_n,dh,hw,qr_temp,dif_p,dif_n)
do k = 1, riv_count

 if(domain_riv_idx(k) .eq. 2) cycle

 zb_p = zb_riv_idx(k)
 hr_p = hr_idx(k)
 dif_p = dif_riv_idx(k)

 distance = dis_riv_idx(k)

 ! information of the destination cell
 kk = down_riv_idx(k)
 zb_n = zb_riv_idx(kk)
 hr_n = hr_idx(kk)
 dif_n = dif_riv_idx(kk)

 ! diffusion wave
 dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance ! diffussion

 ! kinematic wave
 if( dif_p .eq. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

 ! the destination cell is outlet (domain = 2)
 if( domain_riv_idx(kk) .eq. 2 ) dh = (zb_p + hr_p - zb_n) / distance ! kinematic wave (+hr_p)

 ! the destination cell is outlet (domain = 2 ) with water depth boundary water table
 if( domain_riv_idx(kk) .eq. 2 .and. bound_riv_wlev_switch .ge. 1 ) then
  ! ver 1.4.2 mod by T.Sayama on June 24, 2015
  !if( bound_riv_wlev_idx(1, kk) .le. -100.0 ) cycle ! not boundary
  !dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance ! diffussion
  if( bound_riv_wlev_idx(1, kk) .gt. -100.0 ) dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance ! diffussion 
 endif

 ! kinematic wave (for dam and water gate) modified by T.Sayama on Aug 9, 2021
 if ( damflg(k) .gt. 0 ) then
  dh = max((zb_p - zb_n) / distance, 0.001)
 endif
 if ( damflg(kk) .gt. 0 ) then
  if( dam_floodq(damflg(kk)) .gt. 0 ) then
   dh = max((zb_p - zb_n) / distance, 0.001)
  endif
 endif
      
 ! from a tributary (levee height =< 0) to a main river (levee height > 0) : no reverse flow
 !if( height_idx(k) .le. 0.d0 .and. height_idx(kk) .gt. 0.d0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

 if( dh .ge. 0.d0 ) then
  hw = hr_p
  if( zb_p .lt. zb_n ) hw = max(0.d0, zb_p + hr_p - zb_n)
  !call hq_riv(hw, dh, width_idx(k), qr_temp)
  call hq_riv(hw, dh, k, width_idx(k), qr_temp)
  qr_idx(k) = qr_temp
 else
  ! reverse flow
  hw = hr_n
  if( zb_n .lt. zb_p ) hw = max(0.d0, zb_n + hr_n - zb_p)
  dh = abs(dh)
  !call hq_riv(hw, dh, width_idx(k), qr_temp)
  call hq_riv(hw, dh, kk, width_idx(k), qr_temp)
  qr_idx(k) = -qr_temp
 endif

enddo
!$omp end parallel do

end subroutine qr_calc

! water depth and discharge relationship
!subroutine hq_riv(h, dh, w, q)
subroutine hq_riv(h, dh, k, w, q)
use globals
implicit none

real(8) h, dh, w, q
real(8) a, m, r ! add v1.4.2.4
integer k ! add v1.4

a = sqrt(abs(dh)) / ns_river
!m = 5.d0 / 3.d0
r = (w * h) / (w + 2.d0 * h)
! modified v1.4
!q = a * h ** m * w
! modified v1.4.2.4
q = a * r ** (2.d0 / 3.d0) * w * h

if( sec_map_idx(k) .gt. 0 ) then
 call sec_hq_riv(h, dh, k, q)
endif

! discharge per unit area
! (q multiply by width and divide by area)
!q = q * w / area

end subroutine hq_riv
