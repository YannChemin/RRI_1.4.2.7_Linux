! infiltration
subroutine infilt(hs_idx, gampt_ff_idx, gampt_f_idx)
use globals
implicit none

real(8) hs_idx(slo_count), gampt_ff_idx(slo_count), gampt_f_idx(slo_count)
real(8) gampt_ff_temp
integer k

do k = 1, slo_count

 gampt_f_idx(k) = 0.d0
 gampt_ff_temp = gampt_ff_idx(k)
 if( gampt_ff_temp .le. 0.01d0 ) gampt_ff_temp = 0.01d0

 ! gampt_f_idx(k) : infiltration capacity [m/s]
 ! gampt_ff : accumulated infiltration depth [m]
 gampt_f_idx(k) = ksv_idx(k) * (1.d0 + faif_idx(k) * gammaa_idx(k) / gampt_ff_temp)

 ! gampt_f_idx(k) : infiltration capacity -> infiltration rate [m/s]
 if( gampt_f_idx(k) .ge. hs_idx(k) / dt ) gampt_f_idx(k) = hs_idx(k) / dt

 ! gampt_ff should not exceeds a certain level
 if( infilt_limit_idx(k) .ge. 0.d0 .and. gampt_ff_idx(k) .ge. infilt_limit_idx(k) ) gampt_f_idx(k) = 0.d0

 ! update gampt_ff [m]
 gampt_ff_idx(k) = gampt_ff_idx(k) + gampt_f_idx(k) * dt

 ! hs : hs - infiltration rate * dt [m]
 hs_idx(k) = hs_idx(k) - gampt_f_idx(k) * dt
 if( hs_idx(k) .le. 0.d0 ) hs_idx(k) = 0.d0

enddo

end subroutine infilt
