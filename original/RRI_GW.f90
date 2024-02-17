! RRI_GW

! variable definition (slope)

subroutine funcg(hg_idx, fg_idx, qg_idx )
use globals
implicit none

real(8) hg_idx(slo_count), fg_idx(slo_count)
real(8) qg_idx(i4, slo_count)

integer k, l, kk, itemp, jtemp

fg_idx(:) = 0.d0
qg_idx(:,:) = 0.d0

call qg_calc(hg_idx, qg_idx)

! qg_idx > 0 --> discharge flowing out from a cell

!$omp parallel do
do k = 1, slo_count
 if(gammag_idx(k) .gt. 0.d0) fg_idx(k) = (qg_idx(1,k) + qg_idx(2,k) + qg_idx(3,k) + qg_idx(4,k)) / gammag_idx(k)
enddo
!$omp end parallel do

do k = 1, slo_count
 do l = 1, lmax
  if( dif_slo_idx(k) .eq. 0 .and. l .eq. 2 ) exit ! kinematic -> 1-direction
  kk = down_slo_idx(l, k)
  if( dif_slo_idx(k) .eq. 0 ) kk = down_slo_1d_idx(k)
  if( kk .eq. -1 ) cycle
  if(gammag_idx(k) .gt. 0.d0) fg_idx(kk) = fg_idx(kk) - qg_idx(l, k) / gammag_idx(k)
 enddo
enddo

end subroutine funcg


! lateral gw discharge (slope)
subroutine qg_calc(hg_idx, qg_idx)
use globals
implicit none

real(8) hg_idx(slo_count)
real(8) qg_idx(i4, slo_count), qg

integer k, kk, l
real(8) zb_p, hg_p, gammag_p, kg0_p, ksg_p, fpg_p
real(8) zb_n, hg_n, gammag_n, kg0_n, ksg_n, fpg_n
real(8) dh, distance
real(8) len, hw
integer dif_p, dif_n

qg_idx = 0.d0

!$omp parallel do 
!private(kk,zb_p,hg_p,gammag_p,kg0_p,ksg_p,fpg_p,l,distance,len,zb_n,hg_n,gammag_n,kg0_n,ksg_n,fpg_n,dh,dif_p,dif_n,qg)
do k = 1, slo_count

 zb_p = zb_slo_idx(k)
 hg_p = hg_idx(k)
 ksg_p = ksg_idx(k)
 gammag_p = gammag_idx(k)
 kg0_p = kg0_idx(k)
 fpg_p = fpg_idx(k)
 dif_p = dif_slo_idx(k)
 if( ksg_p .le. 0.d0 ) cycle

 ! 8-direction: lmax = 4, 4-direction: lmax = 2
 do l = 1, lmax ! (1: rightC2: down, 3: right down, 4: left down)
  if( dif_p .eq. 0 .and. l .eq. 2 ) exit ! kinematic -> 1-direction
  kk = down_slo_idx(l, k)
  if( dif_p .eq. 0 ) kk = down_slo_1d_idx(k)
  if( kk .eq. -1 ) cycle

  distance = dis_slo_idx(l, k)
  len = len_slo_idx(l, k)
  if( dif_p .eq. 0 ) distance = dis_slo_1d_idx(k)
  if( dif_p .eq. 0 ) len = len_slo_1d_idx(k)

  ! information of the destination cell
  zb_n = zb_slo_idx(kk)
  hg_n = hg_idx(kk)
  gammag_n = gammag_idx(kk)
  kg0_n = kg0_idx(kk)
  ksg_n = ksg_idx(kk)
  fpg_n = fpg_idx(kk)
  dif_n = dif_slo_idx(kk)
  if( ksg_n .le. 0.d0 ) cycle

  ! diffusion wave
  dh = ((zb_p - hg_p) - (zb_n - hg_n)) / distance

  ! 1-direction : kinematic wave
  if( dif_p .eq. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

  ! water coming in or going out?
  if( dh .ge. 0.d0 ) then
   call hg_calc(gammag_p, kg0_p, ksg_p, fpg_p, hg_p, dh, len, qg)
   qg_idx(l,k) = qg
  else
   ! coming in
   dh = abs(dh)
   call hg_calc(gammag_n, kg0_n, ksg_n, fpg_n, hg_n, dh, len, qg)
   qg_idx(l,k) = -qg
  endif

 enddo
enddo
!$omp end parallel do

end subroutine qg_calc


! water depth and gw discharge relationship
subroutine hg_calc(gammag_p, kg0_p, ksg_p, fpg_p, hg_p, dh, len, qg)
use globals
implicit none

real(8) gammag_p, kg0_p, ksg_p, fpg_p, hg_p, dh, len, qg

qg = dh * kg0_p / fpg_p * exp( - fpg_p * hg_p )

! discharge per unit area
! (qg multiply by width and divide by area)
qg = qg * len / area

end subroutine hg_calc


! gw recharge
subroutine gw_recharge(hs_idx, gampt_ff_idx, hg_idx)
use globals
implicit none

real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), hg_idx(slo_count)

real(8) rech, rech_ga, rech_hs
integer k

rech = 0.d0
rech_ga = 0.d0
rech_hs = 0.d0

do k = 1, slo_count

 rech = ksg_idx(k) * dt
 if( hg_idx(k) .lt. 0.d0 ) cycle
 if( hg_idx(k) * gammag_idx(k) .lt. rech ) rech = hg_idx(k) * gammag_idx(k)
 if( rech .lt. gampt_ff_idx(k) .and. ksv_idx(k) .gt. 0.d0 ) then
  rech_ga = rech
  rech_hs = 0.d0
 else
  if( ksv_idx(k) .gt. 0.d0) rech_ga = gampt_ff_idx(k)
  if( rech - rech_ga .lt. hs_idx(k) ) then
   rech_hs = rech - rech_ga
  else
   rech_hs = hs_idx(k)
  endif
 endif
 hs_idx(k) = hs_idx(k) - rech_hs
 if( ksv_idx(k) .gt. 0.d0 ) gampt_ff_idx(k) = gampt_ff_idx(k) - rech_ga
 rech = rech_ga + rech_hs
 if(gammag_idx(k) .gt. 0.d0) hg_idx(k) = hg_idx(k) - rech / gammag_idx(k)
 rech_hs_tsas(k) = rech_hs / dble(dt)

enddo

end subroutine gw_recharge

! gw lose
subroutine gw_lose(hg_idx)
use globals
implicit none

real(8) hg_idx(slo_count)
integer k

do k = 1, slo_count
 if( rgl_idx(k) .gt. 0.d0) hg_idx(k) = hg_idx(k) + rgl_idx(k) / gammag_idx(k) * dt
enddo

end subroutine gw_lose


! gw exfilt
subroutine gw_exfilt(hs_idx, gampt_ff_idx, hg_idx)
use globals
implicit none

real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), hg_idx(slo_count)

real(8) exfilt, exfilt_ga, exfilt_hs
integer k

do k = 1, slo_count

 exfilt = 0.d0
 exfilt_ga = 0.d0
 exfilt_hs = 0.d0

 if( hg_idx(k) .ge. 0.d0 ) cycle
 exfilt = - hg_idx(k) * gammag_idx(k)
 if( infilt_limit_idx(k) .gt. gampt_ff_idx(k) .and. infilt_limit_idx(k) - gampt_ff_idx(k) .ge. exfilt .and. &
    ksv_idx(k) .gt. 0.d0 .and. infilt_limit_idx(k) .gt. 0.d0 ) then
  exfilt_ga = exfilt
  exfilt_hs = 0.d0
 else
  if( ksv_idx(k) .gt. 0.d0 .and. infilt_limit_idx(k) .gt. 0.d0 ) then
   exfilt_ga = infilt_limit_idx(k) - gampt_ff_idx(k)
  endif
  exfilt_hs = exfilt - exfilt_ga
 endif
 hg_idx(k) = 0.d0
 if( ksg_idx(k) .gt. 0.d0 .and. infilt_limit_idx(k) .gt. 0.d0 ) gampt_ff_idx(k) = gampt_ff_idx(k) + exfilt_ga
 hs_idx(k) = hs_idx(k) + exfilt_hs
 exfilt_hs_tsas(k) = exfilt_hs / dble(dt)

enddo

end subroutine gw_exfilt


! initial setting for hg
subroutine hg_init(hg_idx)
use globals
implicit none
real(8) hg_idx(slo_count)

hg_idx(:) = 0.d0

end subroutine hg_init
