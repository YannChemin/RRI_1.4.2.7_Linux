! RRI_DT_Check.f90

! current version of RRI model does not call this subroutine
! this subroutine is only for checking dt characteristics

!******* DT_Check for River ***********************

subroutine dt_check_riv(hr_idx, tt, ddt_chk_riv)
use globals

implicit none

real(8) hr_idx(riv_count), ddt_chk_riv
integer tt

real(8) dt_cfl_idx(riv_count), dt_vnm_idx(riv_count)
real(8) dt_cfl(ny, nx), dt_vnm(ny, nx)

integer k, kk, i, j
real(8) zb_p, hr_p
real(8) zb_n, hr_n
real(8) dh, distance
real(8) hw

real(8) c, alph

character*256 ofile_cfl
character*256 ofile_vnm
character*256 ofile_dt
character*6 t_char

! Output File
call int2char(tt, t_char)
ofile_cfl = "./dt_chk/cfl_riv_" // trim(t_char) // ".out"
ofile_vnm = "./dt_chk/vnm_riv_" // trim(t_char) // ".out"
ofile_dt = "./dt_check_riv.txt"

open(100, file = ofile_cfl)
open(101, file = ofile_vnm)
if(tt.eq.1) open(901, file = ofile_dt)

dt_cfl = 100000.d0
dt_vnm = 100000.d0
dt_cfl_idx = 100000.d0
dt_vnm_idx = 100000.d0

do k = 1, riv_count

 if(domain_riv_idx(k) .eq. 2) cycle

 zb_p = zb_riv_idx(k)
 hr_p = hr_idx(k)

 distance = dis_riv_idx(k)

 ! information of the destination cell
 kk = down_riv_idx(k)
 zb_n = zb_riv_idx(kk)
 hr_n = hr_idx(kk)

 ! diffusion wave
 dh = ((zb_p + hr_p) - (zb_n + hr_n)) / distance ! diffussion
 ! the destination cell is outlet (domain = 2)
 if( domain_riv_idx(kk) .eq. 2 ) dh = (zb_p + hr_p - zb_n) / distance ! kinematic wave

 ! kinematic wave
 !dh = (zb_p - zb_n) / distance ! kinematic wave

 if( dh .ge. 0.d0 ) then
  hw = hr_p
  if( zb_p .lt. zb_n ) hw = max(0.d0, zb_p + hr_p - zb_n)
 else
  ! reverse flow
  hw = hr_n
  if( zb_n .lt. zb_p ) hw = max(0.d0, zb_n + hr_n - zb_p)
  dh = abs(dh)
 endif

 ! CFL Condition
 c = 5.d0 / 3.d0 / ns_river * sqrt(dh) * hw ** (2.d0 / 3.d0)
 dt_cfl_idx(k) = distance / c

 ! von Neuman Condition (when dh < 0.01 => dh = 0.01)
 if( hw .ge. 0.01 ) then
  if( dh .lt. 0.01 ) dh = 0.01 ! dh < 0.01 -> dh = 0.01
  alph = hw ** (5.d0 / 3.d0) / 2.d0 / ns_river / sqrt(dh)
  dt_vnm_idx(k) = distance ** 2.d0 / 4.d0 / alph
 endif

enddo

call sub_riv_idx2ij( dt_cfl_idx, dt_cfl )
call sub_riv_idx2ij( dt_vnm_idx, dt_vnm )

do i = 1, ny
 write(100,'(10000e10.3)') (dt_cfl(i, j), j = 1, nx)
 write(101,'(10000e10.3)') (dt_vnm(i, j), j = 1, nx)
enddo
write(901,'(i7, 3f14.3)') tt, ddt_chk_riv, minval(dt_cfl_idx), minval(dt_vnm_idx)

write(*,*) "dt_min : ", minval(dt_vnm_idx)

close(100)
close(101)

end subroutine dt_check_riv


!******* DT_Check for Slope ***********************

subroutine dt_check_slo(hs_idx, tt, ddt_chk_slo)
use globals

implicit none

real(8) hs_idx(slo_count), ddt_chk_slo
integer tt

real(8) dt_cfl_idx(slo_count), dt_vnm_idx(slo_count)
real(8) dt_cfl(ny, nx), dt_vnm(ny, nx)
real(8) dt_cfl_temp(4), dt_vnm_temp(4)

integer k, kk, l, i, j
real(8) zb_p, hs_p, ns_p, ka_p, da_p, dm_p, b_p
real(8) zb_n, hs_n, ns_n, ka_n, da_n, dm_n, b_n
real(8) dh, distance, width_slo
real(8) lev_p, lev_n
real(8) len, hw

real(8) ns, ka_temp, da_temp, c, alph

character*256 ofile_cfl
character*256 ofile_vnm
character*256 ofile_dt
character*6 t_char

! Output File
call int2char(tt, t_char)
ofile_cfl = "./dt_chk/cfl_slo_" // trim(t_char) // ".out"
ofile_vnm = "./dt_chk/vnm_slo_" // trim(t_char) // ".out"
ofile_dt = "./dt_check_slo.txt"

open(100, file = ofile_cfl)
open(101, file = ofile_vnm)
if(tt.eq.1) open(902, file = ofile_dt)

dt_cfl = 100000.d0
dt_vnm = 100000.d0
dt_cfl_idx = 100000.d0
dt_vnm_idx = 100000.d0

do k = 1, slo_count

 zb_p = zb_slo_idx(k)
 hs_p = hs_idx(k)
 ns_p = ns_slo_idx(k)
 ka_p = ka_idx(k)
 da_p = soildepth_idx(k) * gammaa_idx(k)
 dm_p = soildepth_idx(k) * gammam_idx(k)
 b_p  = beta_idx(k)

 dt_cfl_temp(:) = 10000.d0
 dt_vnm_temp(:) = 10000.d0

 ! 8-direction: lmax = 4, 4-direction: lmax = 2
 do l = 1, lmax ! (1: rightÅC2: down, 3: right down, 4: left down)

  kk = down_slo_idx(l, k)
  if( kk .eq. -1 ) cycle

  distance = dis_slo_idx(l, k)
  len = len_slo_idx(l, k)

  ! information of the destination cell
  zb_n = zb_slo_idx(kk)
  hs_n = hs_idx(kk)
  ns_n = ns_slo_idx(kk)
  ka_n = ka_idx(kk)
  da_n = soildepth_idx(kk) * gammaa_idx(kk)
  dm_n = soildepth_idx(kk) * gammam_idx(kk)
  b_n = beta_idx(kk)

  call h2lev(hs_p, k, lev_p)
  call h2lev(hs_n, kk, lev_n)

  dh = abs((zb_p + lev_p) - (zb_n + lev_n)) / distance

  ! water coming in or going out?
  if( dh .ge. 0.d0 ) then
   ! going out
   hw = hs_p
   if( zb_p .lt. zb_n ) hw = max(0.d0, zb_p + hs_p - zb_n)
   ns = ns_p
   ka_temp = ka_p
   da_temp = da_p
  else
   ! coming in
   hw = hs_n
   dh = abs(dh)
   if( zb_n .lt. zb_p ) hw = max(0.d0, zb_n + hs_n - zb_p)
   ns = ns_n
   ka_temp = ka_n
   da_temp = da_n
  endif

  ! CFL Condition
  if( hw .ge. da_temp ) then
   c = ka_temp * dh + 5.d0 / 3.d0 / ns * sqrt(dh) * ( hw - da_temp ) ** (2.d0 / 3.d0)
  else
   c = ka_temp * dh
  endif
  dt_cfl_temp(l) = distance / c

  ! von Neuman Condition
  if( hw .ge. 0.01 ) then
   if( dh .lt. 0.01 ) dh = 0.01 ! dh < 0.01 -> dh = 0.01
   alph = hw ** (5.d0 / 3.d0) / 2.d0 / ns_river / sqrt(dh)
   dt_vnm_temp(l) = distance ** 2.d0 / 4.d0 / alph
  endif

 enddo

 dt_cfl_idx(k) = minval( dt_cfl_temp )
 dt_vnm_idx(k) = minval( dt_vnm_temp )

enddo

call sub_slo_idx2ij( dt_cfl_idx, dt_cfl )
call sub_slo_idx2ij( dt_vnm_idx, dt_vnm )

do i = 1, ny
 write(100,'(10000e10.3)') (dt_cfl(i, j), j = 1, nx)
 write(101,'(10000e10.3)') (dt_vnm(i, j), j = 1, nx)
enddo
write(902,'(i7, 3f14.3)') tt, ddt_chk_slo, minval(dt_cfl_idx), minval(dt_vnm_idx)

close(100)
close(101)

end subroutine dt_check_slo

