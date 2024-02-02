! Output for T-SAS program by RRI
!
subroutine RRI_TSAS(tt, hs_idx, hr_idx, hg_idx, qs_ave_idx, &
                    qr_ave_idx, qg_ave_idx, qp_t_idx)
use globals
use dam_mod
implicit none

integer tt

character*256 tsasfolder, tsasdata
character*256 ofile_tsas_hs_id, ofile_tsas_hr_id, ofile_tsas_outlet
character*256 ofile_tsas_sto, ofile_tsas_flux, ofile_tsas_rain
character*256 ofile_tsas_hg_param

real(8) hs_idx(slo_count), hr_idx(riv_count), hg_idx(slo_count)
real(8) qs_ave_idx(i4, slo_count), qg_ave_idx(i4, slo_count)
real(8) qr_ave_idx(riv_count), qp_t_idx(slo_count)
real(8) vr_temp ! add v1.4

real(8) distance, zb_p, hs_p, zb_n, hs_n, lev_p, lev_n, dh

integer tsasid, dif_p
integer i, j, k, l, kk
character*6 t_char

parameter( tsasfolder = "./tsas/" )
parameter( tsasdata = "./tsas/data/" )
!parameter( tsasfolder = "F:/project/rok100/" )
!parameter( tsasdata = "F:/project/rok100/data/" )

call int2char(tt, t_char)

! Initial Output
if(tt.eq.0) then

 ofile_tsas_hs_id = trim(tsasfolder) // "hs_id.txt"
 ofile_tsas_hr_id = trim(tsasfolder) // "hr_id.txt"
 ofile_tsas_outlet = trim(tsasfolder) // "outlet.txt"
 ofile_tsas_hg_param = trim(tsasfolder) // "param.txt"

 open( 1101, file = ofile_tsas_hs_id )
 open( 1102, file = ofile_tsas_hr_id )
 open( 1103, file = ofile_tsas_outlet )
 open( 1104, file = ofile_tsas_hg_param )

 allocate( hs_id(ny, nx), hr_id(ny, nx) )
 allocate( hs_id_idx(slo_count), hr_id_idx(riv_count) )

 hs_id(:,:) = -999
 hr_id(:,:) = -999

 tsasid = 0

 do k = 1, slo_count
  tsasid = tsasid + 1
  i = slo_idx2i(k)
  j = slo_idx2j(k)
  hs_id(i, j) = tsasid
  hs_id_idx(k) = tsasid
  if( domain_slo_idx(k) .eq. 2 ) write(1103, *) tsasid
 enddo

 do k = 1, riv_count
  tsasid = tsasid + 1
  i = riv_idx2i(k)
  j = riv_idx2j(k)
  hr_id(i, j) = tsasid
  hr_id_idx(k) = tsasid
  if( domain_riv_idx(k) .eq. 2 ) write(1103, *) tsasid
 enddo

 write(1101, '("ncols", i15)') nx
 write(1101, '("nrows", i15)') ny
 write(1101, '("xllcorner", f22.12)') xllcorner
 write(1101, '("yllcorner", f22.12)') yllcorner
 write(1101, '("cellsize", f20.12)') cellsize
 write(1101, '("NODATA_value", i10)') -999

 write(1102, '("ncols", i15)') nx
 write(1102, '("nrows", i15)') ny
 write(1102, '("xllcorner", f22.12)') xllcorner
 write(1102, '("yllcorner", f22.12)') yllcorner
 write(1102, '("cellsize", f20.12)') cellsize
 write(1102, '("NODATA_value", i10)') -999

 do i = 1, ny
  write(1101, 10001) (hs_id(i, j), j = 1, nx)
  write(1102, 10002) (hr_id(i, j), j = 1, nx)
 enddo

 write(1104, '(e21.10e3)') ns_slope(1)
 write(1104, '(e21.10e3)') soildepth(1)
 write(1104, '(e21.10e3)') gammaa(1)
 write(1104, '(e21.10e3)') ka(1)
 write(1104, '(e21.10e3)') gammam(1)
 write(1104, '(e21.10e3)') beta(1)
 write(1104, '(e21.10e3)') ksg(1)
 write(1104, '(e21.10e3)') gammag(1)
 write(1104, '(e21.10e3)') kg0(1)
 write(1104, '(e21.10e3)') fpg(1)
 write(1104, '(e21.10e3)') rgl(1)
 write(1104, '(e21.10e3)') area

 close(1101)
 close(1102)
 close(1103)
 close(1104)

endif


! Data Output

call int2char(tt, t_char)

ofile_tsas_sto = trim(tsasdata) // "tsas_sto_" // trim(t_char) // ".txt"
ofile_tsas_flux = trim(tsasdata) // "tsas_flux_" // trim(t_char) // ".txt"
ofile_tsas_rain = trim(tsasdata) // "tsas_rain_" // trim(t_char) // ".txt"

open(1103, file = ofile_tsas_sto)
open(1104, file = ofile_tsas_flux)
open(1105, file = ofile_tsas_rain)

! storage

do k = 1, slo_count
 write(1103, 10003) hs_idx(k), hg_idx(k), 1
enddo
do k = 1, riv_count
 call hr2vr(hr_idx(k), k, vr_temp)
 write(1103, 10003) vr_temp, 0.d0, 2
enddo

! flux [m3/s]

! qs
do k = 1, slo_count
 do l = 1, lmax
  dif_p = dif_slo_idx(k)
  if( dif_p .eq. 0 .and. l .eq. 2 ) exit ! kinematic -> 1-direction
  if( dif_p .eq. 0 ) kk = down_slo_1d_idx(k)
  kk = down_slo_idx(l, k)
  if( kk .eq. -1 ) cycle

  ! calc dh/dx
  distance = dis_slo_idx(l, k)
  if( dif_p .eq. 0 ) distance = dis_slo_1d_idx(k)
  zb_p = zb_slo_idx(k)
  hs_p = hs_idx(k)
  zb_n = zb_slo_idx(kk)
  hs_n = hs_idx(kk)
  call h2lev(hs_p, k, lev_p)
  call h2lev(hs_n, kk, lev_n)
  dh = ((zb_p + lev_p) - (zb_n + lev_n)) / distance 
  if( dif_p .eq. 0 ) dh = max( (zb_p - zb_n) / distance, 0.001 )

  write(1104, 10004) hs_id_idx(k), hs_id_idx(kk), & 
   (qs_ave_idx(l, k) * area), dh, (qg_ave_idx(l, k) * area), 1
 enddo
enddo

! qr
do k = 1, riv_count
 if( domain_riv_idx(k) .eq. 2) cycle
 kk = down_riv_idx(k)
 write(1104, 10004) hr_id_idx(k), hr_id_idx(kk), qr_ave_idx(k), 0.d0, 0.d0, 2 ! v1.4
enddo

! qrs
do i = 1, ny
 do j = 1, nx
  if( domain(i,j) .eq. 0 ) cycle
  if( riv(i, j) .eq. 0 ) cycle
  write(1104, 10004) hs_id(i, j), hr_id(i, j), (qrs(i, j) * area), 0.d0, 0.d0, 3
 enddo
enddo


! id, rainfall, aevp, infilt

do k = 1, slo_count
 i = slo_idx2i(k)
 j = slo_idx2j(k)
 write(1105, 10005) hs_id_idx(k), (qp_t_idx(k) * area), (aevp_tsas(k) * area), & 
     (gampt_f(i, j) * area), (rech_hs_tsas(k) * area), (exfilt_hs_tsas(k) * area)
enddo


10001 format(100000i7) ! hs_id
10002 format(100000i7) ! hr_id
10003 format(2e21.10e3, i3) ! sto
10004 format(2i7, 3e21.10e3, i3) ! flux
10005 format(i7, 6e21.10e3) ! rain

close(1103)
close(1104)
close(1105)

end subroutine RRI_TSAS
