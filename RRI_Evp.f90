! RRI_Evp

! Evapotranspiration
subroutine evp( hs_idx, gampt_ff_idx )
use globals
implicit none

real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), qe_t_temp
integer i, j, k, itemp, jtemp

itemp = -1
do jtemp = 1, tt_max_evp
 if( t_evp(jtemp-1) .lt. time .and. time .le. t_evp(jtemp) ) itemp = jtemp
enddo
do i = 1, ny
 do j = 1, nx
  qe_t(i, j) = qe(itemp, evp_i(i), evp_j(j))
 enddo
enddo
call sub_slo_ij2idx( qe_t, qe_t_idx )

! evp_switch = 1 : Allow ET from gampt_ff_idx
!            = 2 : Do not take ET from gampt_ff_idx

do k = 1, slo_count
 i = slo_idx2i(k)
 j = slo_idx2j(k)
 qe_t_temp = qe_t_idx(k)
 
 ! from ver 1.3.3
 !if( dm_idx(k) .gt. 0.d0 .and. hs_idx(k) .lt. dm_idx(k) ) then
 ! qe_t_temp = qe_t_idx(k) * hs_idx(k) / dm_idx(k)
 !endif

 if( evp_switch .eq. 1 ) then
  aevp(i, j) = min( qe_t_temp, (hs_idx(k) + gampt_ff_idx(k)) / dt )
 else
  aevp(i, j) = min( qe_t_temp, hs_idx(k) / dt )
 endif
 aevp_tsas(k) = min( qe_t_temp, hs_idx(k) / dt )

 hs_idx(k) = hs_idx(k) - aevp(i, j) * dt
 if( hs_idx(k) .lt. 0.d0 ) then

  if( evp_switch .eq. 1 ) gampt_ff_idx(k) = gampt_ff_idx(k) + hs_idx(k)

  hs_idx(k) = 0.d0
  if( gampt_ff_idx(k) .lt. 0.d0 ) then
   gampt_ff_idx(k) = 0.d0
  endif
 endif
 aevp_sum = aevp_sum + aevp(i, j) * dt * area
 pevp_sum = pevp_sum + qe_t_idx(k) * dt * area
enddo

end subroutine evp
