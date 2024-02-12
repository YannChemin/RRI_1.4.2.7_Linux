! RRI.f90
!
! coded by T.Sayama
!
! ver 1.4.2.7
!
program RRI
use globals
use runge_mod
use dam_mod, only: dam_switch, dam_vol_temp
use tecout_mod
implicit none

! variable definition

! topographic variable
real(8) nodata
real(8) x1, x2, y1, y2, d1, d2, d3, d4

! rainfall variable
integer, allocatable :: rain_i(:), rain_j(:)
integer tt_max_rain
integer, allocatable :: t_rain(:)
integer nx_rain, ny_rain
real(8), allocatable :: qp(:,:,:), qp_t(:,:), qp_t_idx(:)

! calculation variable
real(8) rho, total_area
real(8) vr_out

real(8), allocatable :: hs(:,:), hr(:,:), hg(:,:), inith(:,:)
real(8), allocatable :: qs_ave(:,:,:), qg_ave(:,:,:), qr_ave(:,:)

real(8), allocatable :: fs(:), hs_idx(:), fr(:), hr_idx(:), fg(:), hg_idx(:)
real(8), allocatable :: qr_idx(:), qr_ave_idx(:), qr_ave_temp_idx(:)
real(8), allocatable :: vr_idx(:)
real(8), allocatable :: qs_idx(:,:), qs_ave_idx(:,:), qs_ave_temp_idx(:,:)
real(8), allocatable :: qg_idx(:,:), qg_ave_idx(:,:), qg_ave_temp_idx(:,:)
real(8), allocatable :: gampt_ff_idx(:), gampt_f_idx(:)

!real(8), allocatable :: rdummy_dim(:)

! other variable
integer i, j, t, k, ios, itemp, jtemp, tt, ii, jj
integer out_next
real(8) out_dt
real(8) rtemp
real(8) ss, sr, si, sg, sinit, sout
integer idummy
real(8) rdummy
real(8) rain_sum
real(8) distance
real(8) ddt_chk_riv, ddt_chk_slo
character*256 ctemp
character*6 t_char
integer div_org_i, div_org_j, div_dest_i, div_dest_j

! calcHydro
character*256 hydro_file, hydro_hr_file
parameter( hydro_file = "hydro.txt" )
parameter( hydro_hr_file = "hydro_hr.txt" )
integer, allocatable :: hydro_i(:), hydro_j(:)
integer maxhydro

! ro
character*256 ofile_ro, outdir1, outdir2
integer kk, l

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 0: FILE NAME AND PARAMETER SETTING
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call RRI_Read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 1: FILE READING 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! max timestep
maxt = lasth * 3600 / dt

! dem file
open( 10, file = demfile, status = "old" )
read(10,*) ctemp, nx
read(10,*) ctemp, ny
read(10,*) ctemp, xllcorner
read(10,*) ctemp, yllcorner
read(10,*) ctemp, cellsize
read(10,*) ctemp, nodata
close(10)

allocate (zs(ny, nx), zb(ny, nx), zb_riv(ny, nx), domain(ny, nx))
call read_gis_real(demfile, zs)

! flow accumulation file
allocate (riv(ny, nx), acc(ny, nx))
call read_gis_int(accfile, acc)

! flow direction file
allocate (dir(ny, nx))
call read_gis_int(dirfile, dir)

! landuse file
allocate( land(ny, nx) )
land = 1
if( land_switch.eq.1 ) then
 call read_gis_int(landfile, land)
endif

! land : 1 ... num_of_landuse
write(*,*) "num_of_landuse : ", num_of_landuse
where( land .le. 0 .or. land .gt. num_of_landuse ) land = num_of_landuse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 2: CALC PREPARATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! dx, dy calc
! d1: south side length
x1 = xllcorner
y1 = yllcorner
x2 = xllcorner + nx * cellsize
y2 = yllcorner
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d1 )

! d2: north side length
x1 = xllcorner
y1 = yllcorner + ny * cellsize
x2 = xllcorner + nx * cellsize
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d2 )

! d3: west side length
x1 = xllcorner
y1 = yllcorner
x2 = xllcorner
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d3 )

! d4: east side length
x1 = xllcorner + nx * cellsize
y1 = yllcorner
x2 = xllcorner + nx * cellsize
y2 = yllcorner + ny * cellsize
if( utm.eq.0 ) call hubeny_sub( x1, y1, x2, y2, d4 )

if( utm.eq.1 ) then
 dx = cellsize
 dy = cellsize
else
 dx = (d1 + d2) / 2.d0 / real(nx)
 dy = (d3 + d4) / 2.d0 / real(ny)
endif
write(*,*) "dx [m] : ", dx, "dy [m] : ", dy

! length and area of each cell
length = sqrt(dx * dy)
area = dx * dy

! river widhth, depth, leavy height, river length, river area ratio
allocate ( width(ny, nx), depth(ny, nx), height(ny, nx), len_riv(ny, nx), area_ratio(ny, nx) )

width = 0.d0
depth = 0.d0
height = 0.d0
len_riv = 0.d0

area_ratio = 0.d0

riv = 0 ! slope cell
if( riv_thresh .gt. 0 ) then
 where(acc .gt. riv_thresh) riv = 1 ! river cell
endif

where(riv.eq.1) width = width_param_c * ( acc * dx * dy * 1.d-6 ) ** width_param_s
where(riv.eq.1) depth = depth_param_c * ( acc * dx * dy * 1.d-6 ) ** depth_param_s
where(riv.eq.1 .and. acc.gt.height_limit_param) height = height_param

! river data is replaced by the information in files
if( rivfile_switch .ge. 1 ) then
 riv = 0
 riv_thresh = 1
 call read_gis_real(widthfile, width)
 where(width .gt. 0) riv = 1 ! river cells (if width >= 0.)
 call read_gis_real(depthfile, depth)
 call read_gis_real(heightfile, height)
 where( height(:,:) .lt. 0.d0 ) height(:,:) = 0.d0
endif 
where(riv.eq.1) len_riv = length

! river cross section is set by section files
allocate( sec_map(ny, nx) )
sec_map = 0
if( sec_switch.eq.1 ) then
 call read_gis_int(sec_map_file, sec_map)
 sec_id_max = maxval( sec_map(:,:) )
 call set_section
endif
where(riv.eq.1) len_riv = length ! added on Dec. 27, 2021

! river length is set by input file
allocate( sec_length(ny, nx) )
sec_length = 0
if( sec_length_switch.eq.1 ) then
 call read_gis_real(sec_length_file, sec_length)
 where(sec_length .gt. 0.d0) len_riv = sec_length
endif

if( rivfile_switch .eq. 2 ) then
 ! levee on both river and slope grid cells : zs is increased with height
 where( height .gt. 0.d0 ) zs = zs + height
else
 ! levee on only slope grid cells : zs is increased with height
 where( height .gt. 0.d0 .and. riv .eq. 0 ) zs = zs + height
endif

!where(riv.eq.1) area_ratio = width / length
!where(riv.eq.1) area_ratio = width / len_riv
where(riv.eq.1) area_ratio = width * len_riv / area ! modified by T.Sayama on Nov. 27, 2021

zb_riv = zs
do i = 1, ny
 do j = 1, nx
  zb(i, j) = zs(i, j) - soildepth(land(i,j))
  if(riv(i, j) .eq. 1) zb_riv(i, j) = zs(i, j) - depth(i, j)
 enddo
enddo

! domain setting

! domain = 0 : outside the domain
! domain = 1 : inside the domain
! domain = 2 : outlet point (where dir(i,j) = 0 or dir(i,j) = -1),
!              and cells located at edges
domain = 0
num_of_cell = 0
do i = 1, ny
 do j = 1, nx
  if( zs(i, j) .gt. -100.d0 ) then
   domain(i, j) = 1
   if( dir(i, j).eq.0 ) domain(i, j) = 2
   if( dir(i, j).eq.-1 ) domain(i, j) = 2
   num_of_cell = num_of_cell + 1
  endif
 enddo
enddo
write(*,*) "num_of_cell : ", num_of_cell
write(*,*) "total area [km2] : ", num_of_cell * area / (10.d0 ** 6.0d0)

write(*,*) "River index setting"
! river index setting
call riv_idx_setting

write(*,*) "Slope index setting"

! slope index setting
call slo_idx_setting

write(*,*) "Reading Dam File"

! reading dam file
call dam_read

write(*,*) "Allocate Initial conditions"
! initial condition
allocate(hs(ny, nx), hr(ny, nx), hg(ny, nx), gampt_ff(ny, nx))
allocate(gampt_f(ny, nx), qrs(ny, nx))

hr = -0.1d0
hs = -0.1d0
hg = -0.1d0
gampt_ff = 0.d0
gampt_f = 0.d0
qrs = 0.d0

!where(riv.eq.1) hr = init_cond_riv
!where(domain.eq.1) hs = init_cond_slo
where(riv.eq.1) hr = 0.d0
where(domain.eq.1) hs = 0.d0
where(domain.eq.2) hs = 0.d0

! for slope cells
! if init_slo_switch = 1 => read from file

if(init_slo_switch .eq. 1) then
 write(*,*) "Init Slope Switch"

 allocate( inith(ny, nx) )
 inith = 0.d0
 open(13, file = initfile_slo, status = "old")
 do i = 1, ny
  read(13,*) (inith(i,j), j = 1, nx)
 enddo
 where(inith .le. 0.d0) inith = 0.d0
 where(domain.eq.1 .and. inith .ge. 0.d0) hs = inith
 deallocate( inith )
 close(13)

endif

! for river cells
! if init_riv_switch = 1 => read from file

if(init_riv_switch .eq. 1) then
 write(*,*) "Init River Switch"

 allocate( inith(ny, nx) )
 inith = 0.d0
 open(13, file = initfile_riv, status = "old")
 do i = 1, ny
  read(13,*) (inith(i,j), j = 1, nx)
 enddo
 where(inith .le. 0.d0) inith = 0.d0
 where(riv.eq.1 .and. inith .ge. 0.d0) hr = inith
 deallocate( inith )
 close(13)

endif

! for slope cells
! if init_gw_switch = 1 => read from file

if(init_gw_switch .eq. 1) then
 write(*,*) "Init GW Switch"

 allocate( inith(ny, nx) )
 inith = 0.d0
 open(13, file = initfile_gw, status = "old")
 do i = 1, ny
  read(13,*) (inith(i,j), j = 1, nx)
 enddo
 where(inith .le. 0.d0) inith = 0.d0
 where(domain.eq.1 .and. inith .ge. 0.d0) hg = inith
 deallocate( inith )
 close(13)

endif

! if init_gampt_ff_switch = 1 => read from file

if(init_gampt_ff_switch .eq. 1) then
 write(*,*) "Init GAMPT FF Switch"

 allocate( inith(ny, nx) )
 inith = 0.d0
 open(13, file = initfile_gampt_ff, status = "old")
 do i = 1, ny
  read(13,*) (inith(i,j), j = 1, nx)
 enddo
 where(inith .le. 0.d0) inith = 0.d0
 where(domain.eq.1) gampt_ff = inith
 deallocate( inith )
 close(13)

endif

write(*,*) "Boundary conditions"
! boundary conditions
call read_bound

! div file
div_id_max = 0
if( div_switch.eq.1 ) then
 open( 20, file = divfile, status = "old" )
 do
  read(20, *, iostat = ios) div_org_i, div_org_j, div_dest_i, div_dest_j
  if(ios .ne. 0) exit
  div_id_max = div_id_max + 1
 enddo
 write(*,*) "div_id_max : ", div_id_max
 allocate( div_org_idx(div_id_max), div_dest_idx(div_id_max), div_rate(div_id_max) )
 rewind(20)

 do k = 1, div_id_max
  read(20, *) div_org_i, div_org_j, div_dest_i, div_dest_j, div_rate(k)
  div_org_idx(k) = riv_ij2idx( div_org_i, div_org_j )
  div_dest_idx(k) = riv_ij2idx( div_dest_i, div_dest_j )
 enddo
 write(*,*) "done: reading div file"
 close(20)
endif

! emb file
!if( emb_switch.eq.1 ) then
! allocate (emb_r(ny, nx), emb_b(ny, nx))
! call read_gis_real(embrfile, emb_r)
! call read_gis_real(embbfile, emb_b)

! call sub_slo_ij2idx(emb_r, emb_r_idx)
! call sub_slo_ij2idx(emb_b, emb_b_idx)
!endif

! hydro file
if( hydro_switch .eq. 1 ) then
 open( 5, file = location_file, status = "old")
 open(1012, file = hydro_file )
 open(1013, file = hydro_hr_file )
 i = 1
 do
  read(5,*,iostat=ios) ctemp, itemp, itemp
  if(ios.ne.0) exit
  i = i + 1
 enddo
 maxhydro = i-1
 rewind(5)
 allocate( hydro_i(maxhydro), hydro_j(maxhydro) )
 do i = 1, maxhydro
  read(5,*) ctemp, hydro_i(i), hydro_j(i)
 enddo
 close(5)
endif

! dynamic allocation
allocate (qs_ave(i4, ny, nx), qr_ave(ny, nx), qg_ave(i4, ny, nx))

allocate (qr_idx(riv_count), qr_ave_idx(riv_count), qr_ave_temp_idx(riv_count), hr_idx(riv_count))
allocate (fr(riv_count), vr_temp(riv_count), hr_err(riv_count), vr_err(riv_count))
allocate (vr_idx(riv_count))
allocate (kr2(riv_count), kr3(riv_count), kr4(riv_count), kr5(riv_count), kr6(riv_count))

allocate (qs_idx(i4, slo_count), qs_ave_idx(i4, slo_count), qs_ave_temp_idx(i4, slo_count), hs_idx(slo_count))
allocate (qp_t_idx(slo_count))
allocate (fs(slo_count), hs_temp(slo_count), hs_err(slo_count))
allocate (ks2(slo_count), ks3(slo_count), ks4(slo_count), ks5(slo_count), ks6(slo_count))

allocate (qg_idx(i4, slo_count), qg_ave_idx(i4, slo_count), qg_ave_temp_idx(i4, slo_count), hg_idx(slo_count))
allocate (fg(slo_count), hg_temp(slo_count), hg_err(slo_count))
allocate (kg2(slo_count), kg3(slo_count), kg4(slo_count), kg5(slo_count), kg6(slo_count))
allocate (gampt_ff_idx(slo_count), gampt_f_idx(slo_count))

allocate (rain_i(ny), rain_j(nx))
allocate (qe_t_idx(slo_count))
allocate (evp_i(ny), evp_j(nx))
allocate (aevp(ny, nx), aevp_tsas(slo_count), exfilt_hs_tsas(slo_count), rech_hs_tsas(slo_count))

! array initialization
qr_ave(:,:) = 0.d0
qr_idx(:) = 0.d0
qr_ave_idx(:) = 0.d0
qr_ave_temp_idx(:) = 0.d0

hr_idx(:) = 0.d0
vr_idx(:) = 0.d0
fr(:) = 0.d0
hr_err(:) = 0.d0
vr_temp(:) = 0.d0
vr_err(:) = 0.d0
kr2(:) = 0.d0
kr3(:) = 0.d0
kr4(:) = 0.d0
kr5(:) = 0.d0
kr6(:) = 0.d0

qs_ave(:,:,:) = 0.d0
qs_idx(:,:) = 0.d0
qs_ave_idx(:,:) = 0.d0
qs_ave_temp_idx(:,:) = 0.d0
hs_idx(:) = 0.d0
qp_t_idx(:) = 0.d0
fs(:) = 0.d0
hs_temp(:) = 0.d0
hs_err(:) = 0.d0
ks2(:) = 0.d0
ks3(:) = 0.d0
ks4(:) = 0.d0
ks5(:) = 0.d0
ks6(:) = 0.d0

qg_ave(:,:,:) = 0.d0
qg_idx(:,:) = 0.d0
qg_ave_idx(:,:) = 0.d0
qg_ave_temp_idx(:,:) = 0.d0
hg_idx(:) = 0.d0
fg(:) = 0.d0
hg_temp(:) = 0.d0
hg_err(:) = 0.d0
gampt_ff_idx(:) = 0.d0
gampt_f_idx(:) = 0.d0

rain_i(:) = 0
rain_j(:) = 0
ksv(:) = 0.d0
faif(:) = 0.d0
aevp(:,:) = 0.d0
hg_idx(:) = 0.d0
evp_i(:) = 0
evp_j(:) = 0
aevp_tsas(:) = 0.d0
exfilt_hs_tsas(:) = 0.d0
rech_hs_tsas(:) = 0.d0

! gw initial setting
if(init_gw_switch .ne. 1) then
 call hg_init( hg_idx )
 call sub_slo_idx2ij( hg_idx, hg )
endif

! initial storage calculation
rain_sum = 0.d0
aevp_sum = 0.d0
pevp_sum = 0.d0
sout = 0.d0
si = 0.d0
sg = 0.d0
open( 1000, file = outfile_storage )
call storage_calc(hs, hr, hg, ss, sr, si, sg)
sinit = ss + sr + si + sg
write(1000, '(1000e15.7)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
  (rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg

! reading rainfall data
open( 11, file = rainfile, status = 'old' )

tt = 0
do
 read(11, *, iostat = ios) t, nx_rain, ny_rain
 do i = 1, ny_rain
  read(11, *, iostat = ios) (rdummy, j = 1, nx_rain)
 enddo
 if( ios.lt.0 ) exit
 tt = tt + 1
enddo
tt_max_rain = tt - 1

allocate( t_rain(0:tt_max_rain), qp(0:tt_max_rain, ny_rain, nx_rain), qp_t(ny, nx) )
rewind(11)

write(*,*) tt_max_rain, nx_rain, ny_rain

qp = 0.d0
qp_t = 0.d0 ! added by T.Sayamaa on Dec 7, 2022 v1.4.2.7
do tt = 0, tt_max_rain
 read(11, *) t_rain(tt), nx_rain, ny_rain
 do i = 1, ny_rain
  read(11, *) (qp(tt, i, j), j = 1, nx_rain)
 enddo
enddo
! unit convert from (mm/h) to (m/s)
qp = qp / 3600.d0 / 1000.d0

do j = 1, nx
 rain_j(j) = int( (xllcorner + (dble(j) - 0.5d0) * cellsize - xllcorner_rain) / cellsize_rain_x ) + 1
enddo
do i = 1, ny
 rain_i(i) = ny_rain - int( (yllcorner + (dble(ny) - dble(i) + 0.5d0) * cellsize - yllcorner_rain) / cellsize_rain_y )
enddo
close(11)

write(*,*) "done: reading rain file"

! reading evp data
if( evp_switch .ne. 0 ) then
 open( 11, file = evpfile, status = 'old' )

 tt = 0
 do
  read(11, *, iostat = ios) t, nx_evp, ny_evp
  do i = 1, ny_evp
   read(11, *, iostat = ios) (rdummy, j = 1, nx_evp)
  enddo
  if( ios.lt.0 ) exit
  tt = tt + 1
 enddo
 tt_max_evp = tt - 1

 allocate( t_evp(0:tt_max_evp), qe(0:tt_max_evp, ny_evp, nx_evp), qe_t(ny, nx) )
 rewind(11)

 write(*,*) "done: reading evp file"

 qe = 0.d0
 do tt = 0, tt_max_evp
  read(11, *) t_evp(tt), nx_evp, ny_evp
  do i = 1, ny_evp
   read(11, *) (qe(tt, i, j), j = 1, nx_evp)
  enddo
 enddo
 ! unit convert from (mm/h) to (m/s)
 qe = qe / 3600.d0 / 1000.d0

 do j = 1, nx
  evp_j(j) = int( (xllcorner + (dble(j) - 0.5d0) * cellsize - xllcorner_evp) / cellsize_evp_x ) + 1
 enddo
 do i = 1, ny
  evp_i(i) = ny_evp - int( (yllcorner + (dble(ny) - dble(i) + 0.5d0) * cellsize - yllcorner_evp) / cellsize_evp_y )
 enddo
 close(11)
endif

! For TSAS Output (Initial Condition)
call sub_slo_ij2idx( hs, hs_idx )
call sub_riv_ij2idx( hr, hr_idx )
!call RRI_TSAS(0, hs_idx, hr_idx, hg_idx, qs_ave_idx, &
!              qr_ave_idx, qg_ave_idx, qp_t_idx)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! STEP 3: CALCULATION 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rain_sum = 0.d0
aevp_sum = 0.d0
sout = 0.d0

! output timestep
out_dt = dble(maxt) / dble(outnum)
out_dt = max(1.d0, out_dt)
out_next = nint(out_dt)
tt = 0

do t = 1, maxt

 if(mod(t, 1).eq.0) write(*,*) t, "/", maxt

 !******* RIVER CALCULATION ******************************
 if( riv_thresh .lt. 0 ) go to 2

 ! from time = (t - 1) * dt to t * dt
 time = (t - 1) * dt ! (current time)
 ! time step is initially set to be "dt_riv"
 ddt = dt_riv
 ddt_chk_riv = dt_riv

 qr_ave = 0.d0
 qr_ave_idx = 0.d0
 if( dam_switch .eq. 1 ) dam_vol_temp(:) = 0.d0

 ! hr -> hr_idx
 ! Memo: riv_ij2idx must be here. 
 ! hr_idx cannot be replaced within the following do loop.
 call sub_riv_ij2idx( hr, hr_idx )

 ! from hr_idx -> vr_idx
 do k = 1, riv_count
  call hr2vr(hr_idx(k), k, vr_idx(k))
 enddo

 do

  ! "time + ddt" should be less than "t * dt"
  if(time + ddt .gt. t * dt ) ddt = t * dt - time

  ! boundary condition for river (water depth boundary)
  if( bound_riv_wlev_switch .ge. 1 ) then
   itemp = -1
   do jtemp = 1, tt_max_bound_riv_wlev
    if( t_bound_riv_wlev(jtemp-1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_riv_wlev(jtemp) ) itemp = jtemp
   enddo
   do k = 1, riv_count
    if( bound_riv_wlev_idx(itemp, k) .le. -100.0 ) cycle ! not boundary
    hr_idx(k) = bound_riv_wlev_idx(itemp, k)
    call hr2vr(hr_idx(k), k, vr_idx(k))
   enddo
  endif

1 continue
  qr_ave_temp_idx(:) = 0.d0

  ! Adaptive Runge-Kutta 
  ! (1)
  call funcr( vr_idx, fr, qr_idx )
  vr_temp = vr_idx + b21 * ddt * fr
  where(vr_temp .lt. 0) vr_temp = 0.d0
  qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

  ! (2)
  call funcr( vr_temp, kr2, qr_idx )
  vr_temp = vr_idx + ddt * (b31 * fr + b32 * kr2)
  where(vr_temp .lt. 0) vr_temp = 0.d0
  qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

  ! (3)
  call funcr( vr_temp, kr3, qr_idx )
  vr_temp = vr_idx + ddt * (b41 * fr + b42 * kr2 + b43 * kr3)
  where(vr_temp .lt. 0) vr_temp = 0.d0
  qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

  ! (4)
  call funcr( vr_temp, kr4, qr_idx )
  vr_temp = vr_idx + ddt * (b51 * fr + b52 * kr2 + b53 * kr3 + b54 * kr4)
  where(vr_temp .lt. 0) vr_temp = 0.d0
  qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

  ! (5)
  call funcr( vr_temp, kr5, qr_idx )
  vr_temp = vr_idx + ddt * (b61 * fr + b62 * kr2 + b63 * kr3 + b64 * kr4 + b65 * kr5)
  where(vr_temp .lt. 0) vr_temp = 0.d0
  qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

  ! (6)
  call funcr( vr_temp, kr6, qr_idx )
  vr_temp = vr_idx + ddt * (c1 * fr + c3 * kr3 + c4 * kr4 + c6 * kr6)
  where(vr_temp .lt. 0) vr_temp = 0.d0
  qr_ave_temp_idx = qr_ave_temp_idx + qr_idx * ddt

  ! (e)
  vr_err = ddt * (dc1 * fr + dc3 * kr3 + dc4 * kr4 + dc5 * kr5 + dc6 * kr6)

  hr_err(:) = vr_err(:) / (area * area_ratio_idx(:))

  ! error evaluation
  where( domain_riv_idx .eq. 0 ) hr_err = 0
  errmax = maxval( hr_err ) / eps
  !if(errmax.gt.1.d0 .and. ddt .ge. ddt_min_riv) then
  if(errmax.gt.1.d0 .and. ddt .gt. ddt_min_riv) then ! modified on Jan 7, 2021
   ! try smaller ddt
   ddt = max( safety * ddt * (errmax ** pshrnk), 0.5d0 * ddt )
   ddt = max( ddt, ddt_min_riv ) ! added on Jan 7, 2021
   ddt_chk_riv = ddt
   write(*,*) "shrink (riv): ", ddt, errmax, maxloc( vr_err )
   if(ddt.eq.0) stop 'stepsize underflow'
   if(dam_switch .eq. 1 ) dam_vol_temp(:) = 0.d0
   go to 1
  else
   ! modified on Jan 7, 2021
   ! ˆÈ‰º‚ð’Ç‰Á‚µ‚ÄAddt = ddt_min_riv‚Ìê‡‚ÍÅI‚Ìvr_temp‚©‚çqr_idx‚ðÄŒvŽZ‚µ‚Ä‚Ý‚é
   if(ddt .eq. ddt_min_riv) then
    call funcr( vr_temp, kr6, qr_idx )
    qr_ave_temp_idx = qr_idx * ddt * 6.d0
   endif
   if(time + ddt .gt. t * dt ) ddt = t * dt - time
   time = time + ddt
   vr_idx = vr_temp
   qr_ave_idx = qr_ave_idx + qr_ave_temp_idx
  endif
  if(time.ge.t * dt) exit ! finish for this timestep
 enddo
 qr_ave_idx = qr_ave_idx / dble(dt) / 6.d0

 do k = 1, riv_count
  call vr2hr(vr_idx(k), k, hr_idx(k))
 enddo

 ! hr_idx -> hr, qr_ave_idx -> qr_ave
 call sub_riv_idx2ij( hr_idx, hr )
 call sub_riv_idx2ij( qr_ave_idx, qr_ave )

 if( dam_switch.eq.1 ) call dam_checkstate(qr_ave)

 !******* SLOPE CALCULATION ******************************
2 continue

 ! from time = (t - 1) * dt to t * dt
 time = (t - 1) * dt  ! (current time)
 ! time step is initially set to be "dt"
 ddt = dt
 ddt_chk_slo = dt

 qs_ave = 0.d0
 qs_ave_idx = 0.d0

 ! hs -> hs_idx
 ! Memo: slo_ij2idx must be here. 
 ! hs_idx cannot be replaced within the following do loop.
 call sub_slo_ij2idx( hs, hs_idx )
 call sub_slo_ij2idx( gampt_ff, gampt_ff_idx ) ! modified by T.Sayama on June 10, 2017

 do

  if(time + ddt .gt. t * dt ) ddt = t * dt - time

  ! rainfall
  itemp = -1
  do jtemp = 1, tt_max_rain
   if( t_rain(jtemp-1) .lt. (time + ddt) .and. (time + ddt) .le. t_rain(jtemp) ) itemp = jtemp
  enddo
  do i = 1, ny
   if(rain_i(i) .lt. 1 .or. rain_i(i) .gt. ny_rain ) cycle
   do j = 1, nx
    if(rain_j(j) .lt. 1 .or. rain_j(j) .gt. nx_rain ) cycle
    qp_t(i, j) = qp(itemp, rain_i(i), rain_j(j))
   enddo
  enddo
  call sub_slo_ij2idx( qp_t, qp_t_idx )

  ! boundary condition for slope (water depth boundary)
  if( bound_slo_wlev_switch .ge. 1 ) then
   itemp = -1
   do jtemp = 1, tt_max_bound_slo_wlev
    if( t_bound_slo_wlev(jtemp-1) .lt. (time + ddt) .and. (time + ddt) .le. t_bound_slo_wlev(jtemp) ) itemp = jtemp
   enddo
   do k = 1, slo_count
    if( bound_slo_wlev_idx(itemp, k) .le. -100.0 ) cycle ! not boundary
    hs_idx(k) = bound_slo_wlev_idx(itemp, k)
   enddo
  endif

3 continue
  qs_ave_temp_idx(:,:) = 0.d0

  ! Adaptive Runge-Kutta 
  ! (1)
  call funcs( hs_idx, qp_t_idx, fs, qs_idx )
  hs_temp = hs_idx + b21 * ddt * fs
  where(hs_temp .lt. 0) hs_temp = 0.d0
  qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt

  ! (2)
  call funcs( hs_temp, qp_t_idx, ks2, qs_idx )
  hs_temp = hs_idx + ddt * (b31 * fs + b32 * ks2)
  where(hs_temp .lt. 0) hs_temp = 0.d0
  qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt

  ! (3)
  call funcs( hs_temp, qp_t_idx, ks3, qs_idx )
  hs_temp = hs_idx + ddt * (b41 * fs + b42 * ks2 + b43 * ks3)
  where(hs_temp .lt. 0) hs_temp = 0.d0
  qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt

  ! (4)
  call funcs( hs_temp, qp_t_idx, ks4, qs_idx )
  hs_temp = hs_idx + ddt * (b51 * fs + b52 * ks2 + b53 * ks3 + b54 * ks4)
  where(hs_temp .lt. 0) hs_temp = 0.d0
  qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt

  ! (5)
  call funcs( hs_temp, qp_t_idx, ks5, qs_idx )
  hs_temp = hs_idx + ddt * (b61 * fs + b62 * ks2 + b63 * ks3 + b64 * ks4 + b65 * ks5)
  where(hs_temp .lt. 0) hs_temp = 0.d0
  qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt

  ! (6)
  call funcs( hs_temp, qp_t_idx, ks6, qs_idx )
  hs_temp = hs_idx + ddt * (c1 * fs + c3 * ks3 + c4 * ks4 + c6 * ks6)
  where(hs_temp .lt. 0) hs_temp = 0.d0
  qs_ave_temp_idx = qs_ave_temp_idx + qs_idx * ddt

  ! (e)
  hs_err = ddt * (dc1 * fs + dc3 * ks3 + dc4 * ks4 + dc5 * ks5 + dc6 * ks6)

  ! error evaluation
  where( domain_slo_idx .eq. 0 ) hs_err = 0.d0
  errmax = maxval( hs_err ) / eps

  !if(errmax.gt.1.d0 .and. ddt .ge. ddt_min_slo) then
  if(errmax.gt.1.d0 .and. ddt .gt. ddt_min_slo) then ! modified on Jan 7, 2021
   ! try smaller ddt
   ddt = max( safety * ddt * (errmax ** pshrnk), 0.5d0 * ddt )
   ddt = max( ddt, ddt_min_slo ) ! added on Jan 7, 2021
   ddt_chk_slo = ddt
   write(*,*) "shrink (slo): ", ddt, errmax, maxloc( hs_err )
   if(ddt.eq.0) stop 'stepsize underflow'
   go to 3
  else
   ! "time + ddt" should be less than "t * dt"
   if(time + ddt .gt. t * dt ) ddt = t * dt - time
   time = time + ddt
   hs_idx = hs_temp
   qs_ave_idx = qs_ave_idx + qs_ave_temp_idx
  endif

  ! cumulative rainfall
  do i = 1, ny
   do j = 1, nx
    if( domain(i,j) .ne. 0 ) rain_sum = rain_sum + dble(qp_t(i, j) * area * ddt)
   enddo
  enddo

  if(time.ge.t * dt) exit ! finish for this timestep
 enddo
 qs_ave_idx = qs_ave_idx / dble(dt) / 6.d0 ! modified on ver 1.4.1

 !******* GW CALCULATION ******************************
 if( gw_switch .eq. 0 ) go to 6

 ! from time = (t - 1) * dt to t * dt
 time = (t - 1) * dt  ! (current time)
 ! time step is initially set to be "dt"
 ddt = dt
 ddt_chk_slo = dt

 qg_ave = 0.d0
 qg_ave_idx = 0.d0

 ! hg -> hg_idx
 ! Memo: slo_ij2idx must be here. 
 ! hg_idx cannot be replaced within the following do loop.
 call sub_slo_ij2idx( hg, hg_idx )
 !call sub_slo_ij2idx( gampt_ff, gampt_ff_idx ) ! modified by T.Sayama on June 10, 2017

 ! GW Recharge
 call gw_recharge( hs_idx, gampt_ff_idx, hg_idx )

 ! GW Lose
 call gw_lose( hg_idx )

 do

  if(time + ddt .gt. t * dt ) ddt = t * dt - time

5 continue
  qg_ave_temp_idx(:,:) = 0.d0

  ! Adaptive Runge-Kutta 
  ! (1)
  call funcg( hg_idx, fg, qg_idx )
  hg_temp = hg_idx + b21 * ddt * fg
  qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

  ! (2)
  call funcg( hg_temp, kg2, qg_idx )
  hg_temp = hg_idx + ddt * (b31 * fg + b32 * kg2)
  qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

  ! (3)
  call funcg( hg_temp, kg3, qg_idx )
  hg_temp = hg_idx + ddt * (b41 * fg + b42 * kg2 + b43 * kg3)
  qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

  ! (4)
  call funcg( hg_temp, kg4, qg_idx )
  hg_temp = hg_idx + ddt * (b51 * fg + b52 * kg2 + b53 * kg3 + b54 * kg4)
  qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

  ! (5)
  call funcg( hg_temp, kg5, qg_idx )
  hg_temp = hg_idx + ddt * (b61 * fg + b62 * kg2 + b63 * kg3 + b64 * kg4 + b65 * kg5)
  qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

  ! (6)
  call funcg( hg_temp, kg6, qg_idx )
  hg_temp = hg_idx + ddt * (c1 * fg + c3 * kg3 + c4 * kg4 + c6 * kg6)
  qg_ave_temp_idx = qg_ave_temp_idx + qg_idx * ddt

  ! (e)
  hg_err = ddt * (dc1 * fg + dc3 * kg3 + dc4 * kg4 + dc5 * kg5 + dc6 * kg6)

  ! error evaluation
  where( domain_slo_idx .eq. 0 ) hg_err = 0.d0
  errmax = maxval( hg_err ) / eps

  !if(errmax.gt.1.d0 .and. ddt .ge. ddt_min_slo) then
  if(errmax.gt.1.d0 .and. ddt .gt. ddt_min_slo) then ! modified on Jan 7, 2021
   ! try smaller ddt
   ddt = max( safety * ddt * (errmax ** pshrnk), 0.5d0 * ddt )
   ddt = max( ddt, ddt_min_slo ) ! added on Jan 7, 2021
   ddt_chk_slo = ddt
   write(*,*) "shrink (gw): ", ddt, errmax, maxloc( hg_err )
   if(ddt.eq.0) stop 'stepsize underflow'
   go to 5
  else
   ! "time + ddt" should be less than "t * dt"
   if(time + ddt .gt. t * dt ) ddt = t * dt - time
   time = time + ddt
   hg_idx = hg_temp
   qg_ave_idx = qg_ave_idx + qg_ave_temp_idx
  endif

  if(time.ge.t * dt) exit ! finish for this timestep
 enddo
 qg_ave_idx = qg_ave_idx / dble(dt) / 6.d0

 time = t * dt

 !******* GW Exfiltration ********************************
 call gw_exfilt( hs_idx, gampt_ff_idx, hg_idx )

6 continue

 !******* Evapotranspiration *****************************
 if( evp_switch .ne. 0 ) call evp(hs_idx, gampt_ff_idx)

 ! hs_idx -> hs
 call sub_slo_idx2ij( hs_idx, hs )
 call sub_slo_idx2ij4( qs_ave_idx, qs_ave )
 call sub_slo_idx2ij( hg_idx, hg )
 call sub_slo_idx2ij4( qg_ave_idx, qg_ave )
 call sub_slo_idx2ij( gampt_ff_idx, gampt_ff )

 !******* LEVEE BREAK ************************************
 !call levee_break(t, hr, hs, xllcorner, yllcorner, cellsize)

 !******* RIVER-SLOPE INTERACTIONS ***********************
 if( riv_thresh .ge. 0 ) call funcrs(hr, hs)
 call sub_riv_ij2idx( hr, hr_idx )
 call sub_slo_ij2idx( hs, hs_idx )

 !******* INFILTRATION (Green Ampt) **********************
 call infilt(hs_idx, gampt_ff_idx, gampt_f_idx)
 call sub_slo_idx2ij( hs_idx, hs )
 call sub_slo_idx2ij( gampt_ff_idx, gampt_ff )
 call sub_slo_idx2ij( gampt_f_idx, gampt_f )

 !******* SET WATER DEPTH 0 AT DOMAIN = 2 ****************
 do i = 1, ny
  do j = 1, nx
   if( domain(i,j) .eq. 2  ) then
    sout = sout + hs(i,j) * area
    hs(i,j) = 0.d0
    if( riv(i,j) .eq. 1) then
     call hr2vr(hr(i, j), riv_ij2idx(i,j), vr_out)
     sout = sout + vr_out
     hr(i,j) = 0.d0
    endif
   endif
  enddo
 enddo

 ! hs -> hs_idx, hr -> hr_idx, hg -> hg_idx
 call sub_riv_ij2idx( hr, hr_idx )
 call sub_slo_ij2idx( hs, hs_idx )
 call sub_slo_ij2idx( hg, hg_idx )

 write(*,*) "max hr: ", maxval(hr), "loc : ", maxloc(hr)
 write(*,*) "max hs: ", maxval(hs), "loc : ", maxloc(hs)
 if(gw_switch .eq. 1) write(*,*) "max hg: ", maxval(hg), "loc : ", maxloc(hg)

 !******* OUTPUT *****************************************

 ! For TSAS Output
 !call RRI_TSAS(t, hs_idx, hr_idx, hg_idx, qs_ave_idx, &
 !              qr_ave_idx, qg_ave_idx, qp_t_idx)

 if( hydro_switch .eq. 1 .and. mod(int(time), 3600) .eq. 0 ) write(1012, '(f12.2, 10000f14.5)') time, &
  (qr_ave(hydro_i(k), hydro_j(k)), k = 1, maxhydro)

 if( hydro_switch .eq. 1 .and. mod(int(time), 3600) .eq. 0 ) write(1013, '(f12.2, 10000f14.5)') time, &
  (hr(hydro_i(k), hydro_j(k)), k = 1, maxhydro) ! added by T.Sayama on July 1, 2021

 ! open output files
 if( t .eq. out_next ) then

  write(*,*) "OUTPUT :", t, time

  tt = tt + 1
  out_next = nint((tt+1) * out_dt)
  call int2char(tt, t_char)

  where(domain .eq. 0) hs = -0.1d0
  if(riv_thresh.ge.0) where(domain .eq. 0) hr = -0.1d0
  if(riv_thresh.ge.0) where(domain .eq. 0) qr_ave = -0.1d0
  where(domain .eq. 0) gampt_ff = -0.1d0
  where(domain .eq. 0) aevp = -0.1d0
  if( evp_switch .ne. 0 ) where(domain .eq. 0) qe_t = -0.1d0
  where(domain .eq. 0) hg = -0.1d0

  if(outswitch_hs .eq. 1) ofile_hs = trim(outfile_hs) // trim(t_char) // ".out"
  if(outswitch_hs .eq. 2) ofile_hs = trim(outfile_hs) // trim(t_char) // ".bin"
  if(outswitch_hr .eq. 1) ofile_hr = trim(outfile_hr) // trim(t_char) // ".out"
  if(outswitch_hr .eq. 2) ofile_hr = trim(outfile_hr) // trim(t_char) // ".bin"
  if(outswitch_hg .eq. 1) ofile_hg = trim(outfile_hg) // trim(t_char) // ".out"
  if(outswitch_hg .eq. 2) ofile_hg = trim(outfile_hg) // trim(t_char) // ".bin"
  if(outswitch_qr .eq. 1) ofile_qr = trim(outfile_qr) // trim(t_char) // ".out"
  if(outswitch_qr .eq. 2) ofile_qr = trim(outfile_qr) // trim(t_char) // ".bin"
  if(outswitch_qu .eq. 1) ofile_qu = trim(outfile_qu) // trim(t_char) // ".out"
  if(outswitch_qu .eq. 2) ofile_qu = trim(outfile_qu) // trim(t_char) // ".bin"
  if(outswitch_qv .eq. 1) ofile_qv = trim(outfile_qv) // trim(t_char) // ".out"
  if(outswitch_qv .eq. 2) ofile_qv = trim(outfile_qv) // trim(t_char) // ".bin"
  if(outswitch_gu .eq. 1) ofile_gu = trim(outfile_gu) // trim(t_char) // ".out"
  if(outswitch_gu .eq. 2) ofile_gu = trim(outfile_gu) // trim(t_char) // ".bin"
  if(outswitch_gv .eq. 1) ofile_gv = trim(outfile_gv) // trim(t_char) // ".out"
  if(outswitch_gv .eq. 2) ofile_gv = trim(outfile_gv) // trim(t_char) // ".bin"
  if(outswitch_gampt_ff .eq. 1) ofile_gampt_ff = trim(outfile_gampt_ff) // trim(t_char) // ".out"
  if(outswitch_gampt_ff .eq. 2) ofile_gampt_ff = trim(outfile_gampt_ff) // trim(t_char) // ".bin"

  if(outswitch_hs .eq. 1) open( 100, file = ofile_hs )
  if(outswitch_hr .eq. 1) open( 101, file = ofile_hr )
  if(outswitch_hg .eq. 1) open( 102, file = ofile_hg )
  if(outswitch_qr .eq. 1) open( 103, file = ofile_qr )
  if(outswitch_qu .eq. 1) open( 104, file = ofile_qu )
  if(outswitch_qv .eq. 1) open( 105, file = ofile_qv )
  if(outswitch_gu .eq. 1) open( 106, file = ofile_gu )
  if(outswitch_gv .eq. 1) open( 107, file = ofile_gv )
  if(outswitch_gampt_ff .eq. 1) open( 108, file = ofile_gampt_ff )

  if(outswitch_hs .eq. 2) open( 100, file = ofile_hs, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_hr .eq. 2) open( 101, file = ofile_hr, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_hg .eq. 2) open( 102, file = ofile_hr, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_qr .eq. 2) open( 103, file = ofile_qr, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_qu .eq. 2) open( 104, file = ofile_qu, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_qv .eq. 2) open( 105, file = ofile_qv, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_gu .eq. 2) open( 106, file = ofile_gu, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_gv .eq. 2) open( 107, file = ofile_gv, form = 'unformatted', access = 'direct', recl = nx*ny*4 )
  if(outswitch_gampt_ff .eq. 2) open( 108, file = ofile_gampt_ff, form = 'unformatted', access = 'direct', recl = nx*ny*4 )

  ! output (ascii)
  do i = 1, ny
   if(outswitch_hs .eq. 1) write(100,'(10000f14.5)') (hs(i, j), j = 1, nx)
   if(outswitch_hr .eq. 1) write(101,'(10000f14.5)') (hr(i, j), j = 1, nx)
   !if(outswitch_hr .eq. 1) write(101,'(10000f14.5)') ((hr(i, j) + zb_riv(i, j)), j = 1, nx)
   if(outswitch_hg .eq. 1) write(102,'(10000f14.5)') (hg(i, j), j = 1, nx)
   if(outswitch_qr .eq. 1) write(103,'(10000f14.5)') ((qr_ave(i,j)), j = 1, nx) ! [m3/s]
   if(outswitch_qu .eq. 1) write(104,'(10000f14.5)') &
     (((qs_ave(1, i, j) + (qs_ave(3, i, j) - qs_ave(4, i, j)) / 2.d0) * area), j = 1, nx)
   !if(outswitch_qv .eq. 1) write(105,'(10000f14.5)') &
   ! (((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j)) / 2.d0) * area), j = 1, nx)
   if(outswitch_qv .eq. 1) write(105,'(10000e14.5)') &
     (((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j)) / 2.d0) * area), j = 1, nx)
   if(outswitch_gu .eq. 1) write(106,'(10000f14.8)') &
     (((qg_ave(1, i, j) + (qg_ave(3, i, j) - qg_ave(4, i, j)) / 2.d0) * area), j = 1, nx)
   if(outswitch_gv .eq. 1) write(107,'(10000f14.5)') &
     (((qg_ave(2, i, j) + (qg_ave(3, i, j) + qg_ave(4, i, j)) / 2.d0) * area), j = 1, nx)
   if(outswitch_gampt_ff .eq. 1) write(108,'(10000f14.5)') (gampt_ff(i, j), j = 1, nx)
  enddo

  ! output (binary)
  if(outswitch_hs .eq. 2) write(100,rec=1) ((hs(i,j), j = 1, nx), i = ny, 1, -1)
  if(outswitch_hr .eq. 2) write(101,rec=1) ((hr(i,j), j = 1, nx), i = ny, 1, -1)
  if(outswitch_hg .eq. 2) write(102,rec=1) ((hg(i,j), j = 1, nx), i = ny, 1, -1)
  if(outswitch_qr .eq. 2) write(103,rec=1) ((qr_ave(i,j), j = 1, nx), i = ny, 1, -1) ! [m3/s]
  if(outswitch_qu .eq. 2) write(104,rec=1) (((qs_ave(1, i, j) + (qs_ave(3, i, j) - qs_ave(4, i, j)) / 2.d0) * area), &
i = ny, 1, -1)
  if(outswitch_qv .eq. 2) write(105,rec=1) (((qs_ave(2, i, j) + (qs_ave(3, i, j) + qs_ave(4, i, j)) / 2.d0) * area), &
i = ny, 1, -1)
  if(outswitch_gu .eq. 2) write(106,rec=1) (((qg_ave(1, i, j) + (qg_ave(3, i, j) - qg_ave(4, i, j)) / 2.d0) * area), &
i = ny, 1, -1)
  if(outswitch_gv .eq. 2) write(107,rec=1) (((qg_ave(2, i, j) + (qg_ave(3, i, j) + qg_ave(4, i, j)) / 2.d0) * area), &
i = ny, 1, -1)
  if(outswitch_gampt_ff .eq. 2) write(108,rec=1) ((gampt_ff(i,j), j = 1, nx), i = ny, 1, -1)

  if(outswitch_hs .ne. 0) close(100)
  if(outswitch_hr .ne. 0) close(101)
  if(outswitch_hg .ne. 0) close(102)
  if(outswitch_qr .ne. 0) close(103)
  if(outswitch_qu .ne. 0) close(104)
  if(outswitch_qv .ne. 0) close(105)
  if(outswitch_gu .ne. 0) close(106)
  if(outswitch_gv .ne. 0) close(107)
  if(outswitch_gampt_ff .ne. 0) close(108)

  if( tec_switch .eq. 1 ) then
   if (tt .eq. 1) then
    call Tecout_alloc(nx, ny, 4)
    call Tecout_mkGrid(dx, dy, zs)
    call Tecout_write_initialize(tt, width, depth, height, area_ratio)
   endif
   call Tecout_write(tt, qp_t, hr, qr_ave, hs, area)
  end if

  ! For dt_check
  !call dt_check_riv(hr_idx, tt, ddt_chk_riv)
  !call dt_check_slo(hs_idx, tt, ddt_chk_slo)

  if( dam_switch .eq. 1 ) then
   if( tt .eq. 1 ) open(1001, file = "./out/dam_out.txt")
   call dam_write
   if(t .eq. maxt ) then
    open(1002, file = "./out/damcnt_out.txt")
    call dam_write_cnt
   endif
  endif

  !******* OUTPUT FOR UNSTEADY MODEL (DIR = -1) *********** added on Sep 15, 2019
  itemp = 0
  do k = 1, riv_count
   i = riv_idx2i(k)
   j = riv_idx2j(k)
   if(dir(i, j) .eq. -1) itemp = 1
  enddo

  if( itemp .eq. 1) then
   ofile_ro = './out/ro_' // trim(t_char) // ".out"
   open(111, file = ofile_ro)
   do k = 1, riv_count
    i = riv_idx2i(k)
    j = riv_idx2j(k)
    kk = down_riv_idx(k)
    if(k.eq.kk) cycle
    ii = riv_idx2i(kk)
    jj = riv_idx2j(kk)
    if(dir(ii,jj) .eq. -1) then
     write(111,'(5i6, 10000f14.5)') 1, i, j, ii, jj, qr_ave_idx(k)
    endif
   enddo

   do k = 1, slo_count
    i = slo_idx2i(k)
    j = slo_idx2j(k)
    if(dir(i, j).eq.-1) cycle
    do l = 1, lmax
     kk = down_slo_idx(l, k)
     if(kk.le.0) cycle
     ii = slo_idx2i(kk)
     jj = slo_idx2j(kk)
     if(dir(ii,jj) .eq. -1) then
      write(111,'(5i6, 10000f14.5)') 2, i, j, ii, jj, qs_ave(l, i, j) * area
     endif
    enddo
   enddo

   do k = 1, slo_count
    i = slo_idx2i(k)
    j = slo_idx2j(k)
    if(dir(i,j) .eq. -1) then
     do l = 1, lmax
      kk = down_slo_idx(l, k)
      if(kk.le.0) cycle
      ii = slo_idx2i(kk)
      jj = slo_idx2j(kk)
      if(domain(ii, jj) .ne. 1) cycle
      write(111,'(5i6, 10000f14.5)') 3, ii, jj, i, j, -1.d0 * qs_ave(l, ii, jj) * area
     enddo
    endif
   enddo

   ! added on Feb 14, 2021 to include rainfall on grid cells with dir = -1
   do k = 1, slo_count
    i = slo_idx2i(k)
    j = slo_idx2j(k)
    if(dir(i,j) .eq. -1) then
     write(111,'(5i6, 10000f14.5)') 4, i, j, i, j, qp_t_idx(k) * area
    endif
   enddo

   close(111)
  endif
 endif

 ! check water balance
 if(mod(t, 1).eq.0) then
  call storage_calc(hs, hr, hg, ss, sr, si, sg)
  write(*, '(6e15.3)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
(rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit)
  write(1000, '(1000e15.7)') rain_sum, pevp_sum, aevp_sum, sout, ss + sr + si + sg, &
(rain_sum - aevp_sum - sout - (ss + sr + si + sg) + sinit), ss, sr, si, sg
 endif

enddo

!pause

end program RRI
