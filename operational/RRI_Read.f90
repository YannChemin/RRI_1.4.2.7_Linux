subroutine RRI_Read

use globals
use dam_mod
use tecout_mod
implicit none

integer i
character*256 format_version

open(1, file = "RRI_Input.txt", status = 'old')

read(1,'(a)') format_version
!write(*,'("format_version : ", a)') trim(adjustl(format_version))

if( format_version .ne. "RRI_Input_Format_Ver1_4_2" ) stop "This RRI model requires RRI_Input_Format_Ver1_4_2"

read(1,*)
!write(*,*)

read(1,'(a)') rainfile
read(1,'(a)') demfile
read(1,'(a)') accfile
read(1,'(a)') dirfile

!write(*,'("rainfile : ", a)') trim(adjustl(rainfile))
!write(*,'("demfile : ", a)') trim(adjustl(demfile))
!write(*,'("accfile : ", a)') trim(adjustl(accfile))
!write(*,'("dirfile : ", a)') trim(adjustl(dirfile))

read(1,*)
!write(*,*)

read(1,*) utm
read(1,*) eight_dir
read(1,*) lasth
read(1,*) dt
read(1,*) dt_riv
read(1,*) outnum
read(1,*) xllcorner_rain
read(1,*) yllcorner_rain
read(1,*) cellsize_rain_x, cellsize_rain_y

!write(*,'("utm : ", i5)') utm
!write(*,'("eight_dir : ", i5)') eight_dir
!write(*,'("lasth : ", i8)') lasth
!write(*,'("dt : ", i12)') dt
!write(*,'("dt_riv : ", i8)') dt_riv
!write(*,'("outnum : ", i8)') outnum
!write(*,'("xllcorner_rain : ", f15.5)') xllcorner_rain
!write(*,'("yllcorner_rain : ", f15.5)') yllcorner_rain
!write(*,'("cellsize_rain_x : ", f15.5, "  cellsize_rain_y : ", f15.5)') cellsize_rain_x, cellsize_rain_y

read(1,*)
!write(*,*)

read(1,*) ns_river
read(1,*) num_of_landuse

allocate( dif(num_of_landuse) )
allocate( ns_slope(num_of_landuse), soildepth(num_of_landuse) )
allocate( gammaa(num_of_landuse) )

read(1,*) (dif(i), i = 1, num_of_landuse)
read(1,*) (ns_slope(i), i = 1, num_of_landuse)
read(1,*) (soildepth(i), i = 1, num_of_landuse)
read(1,*) (gammaa(i), i = 1, num_of_landuse)

!write(*,'("ns_river : ", f12.3)') ns_river
!write(*,'("num_of_landuse : ", i5)') num_of_landuse
!write(*,'("dif : ", 100i5)') (dif(i), i = 1, num_of_landuse)
!write(*,'("ns_slope : ", 100f12.3)') (ns_slope(i), i = 1, num_of_landuse)
!write(*,'("soildepth : ", 100f12.3)') (soildepth(i), i = 1, num_of_landuse)
!write(*,'("gammaa : ", 100f12.3)') (gammaa(i), i = 1, num_of_landuse)

read(1,*)
!write(*,*)

allocate( ksv(num_of_landuse), faif(num_of_landuse) )

read(1,*) (ksv(i), i = 1, num_of_landuse)
read(1,*) (faif(i), i = 1, num_of_landuse)

!write(*,'("ksv : ", 100e12.3)') (ksv(i), i = 1, num_of_landuse)
!write(*,'("faif : ", 100f12.3)') (faif(i), i = 1, num_of_landuse)

read(1,*)
!write(*,*)

allocate( ka(num_of_landuse), gammam(num_of_landuse), beta(num_of_landuse) )

read(1,*) (ka(i), i = 1, num_of_landuse)
read(1,*) (gammam(i), i = 1, num_of_landuse)
read(1,*) (beta(i), i = 1, num_of_landuse)

!write(*,'("ka : ", 100e12.3)') (ka(i), i = 1, num_of_landuse)
!write(*,'("gammam : ", 100f12.3)') (gammam(i), i = 1, num_of_landuse)
!write(*,'("beta : ", 100f12.3)') (beta(i), i = 1, num_of_landuse)

do i = 1, num_of_landuse
 if( gammam(i) .gt. gammaa(i) ) stop "gammam must be smaller than gammaa"
enddo

read(1,*)
!write(*,*)

allocate( ksg(num_of_landuse), gammag(num_of_landuse), kg0(num_of_landuse), fpg(num_of_landuse), rgl(num_of_landuse) )

read(1,*) (ksg(i), i = 1, num_of_landuse)
read(1,*) (gammag(i), i = 1, num_of_landuse)
read(1,*) (kg0(i), i = 1, num_of_landuse)
read(1,*) (fpg(i), i = 1, num_of_landuse)
read(1,*) (rgl(i), i = 1, num_of_landuse)

!write(*,'("ksg : ", 100e12.3)') (ksg(i), i = 1, num_of_landuse)
!write(*,'("gammag : ", 100f12.3)') (gammag(i), i = 1, num_of_landuse)
!write(*,'("kg0 : ", 100e12.3)') (kg0(i), i = 1, num_of_landuse)
!write(*,'("fpg : ", 100f12.3)') (fpg(i), i = 1, num_of_landuse)
!write(*,'("rgl : ", 100e12.3)') (rgl(i), i = 1, num_of_landuse)

read(1,*)
!write(*,*)

read(1,*) riv_thresh
read(1,*) width_param_c
read(1,*) width_param_s
read(1,*) depth_param_c
read(1,*) depth_param_s
read(1,*) height_param
read(1,*) height_limit_param

read(1,*)
!write(*,*)

read(1,*) rivfile_switch
read(1,'(a)') widthfile
read(1,'(a)') depthfile
read(1,'(a)') heightfile

!if(rivfile_switch.eq.0) then
! write(*,'("riv_thresh : ", i7)') riv_thresh
! write(*,'("width_param_c : ", f12.2)') width_param_c
! write(*,'("width_param_s : ", f12.2)') width_param_s
! write(*,'("depth_param_c : ", f12.2)') depth_param_c
! write(*,'("depth_param_s : ", f12.2)') depth_param_s
! write(*,'("height_param : ", f12.2)') height_param
! write(*,'("height_limit_param : ", i10)') height_limit_param
!else
! write(*,'("widthfile : ", a)') trim(adjustl(widthfile))
! write(*,'("depthfile : ", a)') trim(adjustl(depthfile))
! write(*,'("heightfile : ", a)') trim(adjustl(heightfile))
!endif

read(1,*)
!write(*,*)

read(1,*) init_slo_switch, init_riv_switch, init_gw_switch, init_gampt_ff_switch
read(1,"(a)") initfile_slo
read(1,'(a)') initfile_riv
read(1,'(a)') initfile_gw
read(1,'(a)') initfile_gampt_ff

!if(init_slo_switch.ne.0) write(*,'("initfile_slo : ", a)') trim(adjustl(initfile_slo))
!if(init_riv_switch.ne.0) write(*,'("initfile_riv : ", a)') trim(adjustl(initfile_riv))
!if(init_gw_switch.ne.0) write(*,'("initfile_gw : ", a)') trim(adjustl(initfile_gw))
!if(init_gampt_ff_switch.ne.0) write(*,'("initfile_gampt_ff : ", a)') trim(adjustl(initfile_gampt_ff))

read(1,*)
!write(*,*)

read(1,*) bound_slo_wlev_switch, bound_riv_wlev_switch
read(1,'(a)') boundfile_slo_wlev
read(1,'(a)') boundfile_riv_wlev

!if(bound_slo_wlev_switch.ne.0) write(*,'("boundfile_slo_wlev : ", a)') trim(adjustl(boundfile_slo_wlev))
!if(bound_riv_wlev_switch.ne.0) write(*,'("boundfile_riv_wlev : ", a)') trim(adjustl(boundfile_riv_wlev))

read(1,*)
!write(*,*)

read(1,*) bound_slo_disc_switch, bound_riv_disc_switch
read(1,'(a)') boundfile_slo_disc
read(1,'(a)') boundfile_riv_disc

if(bound_slo_disc_switch.ne.0) write(*,'("boundfile_slo_disc : ", a)') trim(adjustl(boundfile_slo_disc))
if(bound_riv_disc_switch.ne.0) write(*,'("boundfile_riv_disc : ", a)') trim(adjustl(boundfile_riv_disc))

read(1,*)
!write(*,*)

read(1,*) land_switch
read(1,'(a)') landfile
if(land_switch.eq.1) write(*,'("landfile : ", a)') trim(adjustl(landfile))

read(1,*)
!write(*,*)

read(1,*) dam_switch
read(1,'(a)') damfile
if(dam_switch.eq.1) write(*,'("damfile : ", a)') trim(adjustl(damfile))

read(1,*)
!write(*,*)

read(1,*) div_switch
read(1,'(a)') divfile
if(div_switch.eq.1) write(*,'("divfile : ", a)') trim(adjustl(divfile))

read(1,*)
!write(*,*)

read(1,*) evp_switch
read(1,'(a)') evpfile
read(1,*) xllcorner_evp
read(1,*) yllcorner_evp
read(1,*) cellsize_evp_x, cellsize_evp_y

!if( evp_switch .ne. 0 ) then
 !write(*,'("evpfile : ", a)') trim(adjustl(evpfile))
 !write(*,'("xllcorner_evp : ", f15.5)') xllcorner_evp
 !write(*,'("yllcorner_evp : ", f15.5)') yllcorner_evp
 !write(*,'("cellsize_evp_x : ", f15.5, " cellsize_evp_y : ", f15.5)') cellsize_evp_x, cellsize_evp_y
!endif

read(1,*)
!write(*,*)

read(1,*) sec_length_switch
read(1,'(a)') sec_length_file
!if(sec_length_switch.eq.1) write(*,'("sec_length : ", a)') trim(adjustl(sec_length_file))

read(1,*)
!write(*,*)

read(1,*) sec_switch
read(1,'(a)') sec_map_file
read(1,'(a)') sec_file
!if(sec_switch.eq.1) write(*,'("sec_map_file : ", a)') trim(adjustl(sec_map_file))
!if(sec_switch.eq.1) write(*,'("sec_file : ", a)') trim(adjustl(sec_file))

read(1,*)
!write(*,*)

!read(1,*) emb_switch
!read(1,'(a)') embrfile
!read(1,'(a)') embbfile
!if(emb_switch.eq.1) write(*,'("embrfile : ", a)') trim(adjustl(embrfile))
!if(emb_switch.eq.1) write(*,'("embbfile : ", a)') trim(adjustl(embbfile))

read(1,*) outswitch_hs, outswitch_hr, outswitch_hg, outswitch_qr, outswitch_qu, outswitch_qv, &
          outswitch_gu, outswitch_gv, outswitch_gampt_ff, outswitch_storage
read(1,'(a)') outfile_hs
read(1,'(a)') outfile_hr
read(1,'(a)') outfile_hg
read(1,'(a)') outfile_qr
read(1,'(a)') outfile_qu
read(1,'(a)') outfile_qv
read(1,'(a)') outfile_gu
read(1,'(a)') outfile_gv
read(1,'(a)') outfile_gampt_ff
read(1,'(a)') outfile_storage

!if(outswitch_hs .ne. 0) write(*,'("outfile_hs : ", a)') trim(adjustl(outfile_hs))
!if(outswitch_hr .ne. 0) write(*,'("outfile_hr : ", a)') trim(adjustl(outfile_hr))
!if(outswitch_hg .ne. 0) write(*,'("outfile_hg : ", a)') trim(adjustl(outfile_hg))
!if(outswitch_qr .ne. 0) write(*,'("outfile_qr : ", a)') trim(adjustl(outfile_qr))
!if(outswitch_qu .ne. 0) write(*,'("outfile_qu : ", a)') trim(adjustl(outfile_qu))
!if(outswitch_qv .ne. 0) write(*,'("outfile_qv : ", a)') trim(adjustl(outfile_qv))
!if(outswitch_gu .ne. 0) write(*,'("outfile_gu : ", a)') trim(adjustl(outfile_gu))
!if(outswitch_gv .ne. 0) write(*,'("outfile_gv : ", a)') trim(adjustl(outfile_gv))
!if(outswitch_gampt_ff .ne. 0) write(*,'("outfile_gampt_ff : ", a)') trim(adjustl(outfile_gampt_ff))
!if(outswitch_storage .ne. 0) write(*,'("outfile_storage : ", a)') trim(adjustl(outfile_storage))

read(1,*)
!write(*,*)

read(1,*) hydro_switch
read(1,'(a)') location_file
!if(hydro_switch .eq. 1) write(*,'("location_file : ", a)') trim(adjustl(location_file))

!write(*,*)

close(1)

! Parameter Check
do i = 1, num_of_landuse
 if( ksv(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 ) &
  stop "Error: both ksv and ka are non-zero."
 if( gammam(i) .gt. gammaa(i) ) &
  stop "Error: gammam must be smaller than gammaa."
enddo

! Set da, dm and infilt_limit
allocate( da(num_of_landuse), dm(num_of_landuse), infilt_limit(num_of_landuse))
da(:) = 0.d0
dm(:) = 0.d0
infilt_limit(:) = 0.d0
do i = 1, num_of_landuse
 if( soildepth(i) .gt. 0.d0 .and. ksv(i) .gt. 0.d0 ) infilt_limit(i) = soildepth(i) * gammaa(i)
 if( soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 ) da(i)= soildepth(i) * gammaa(i)
 if( soildepth(i) .gt. 0.d0 .and. ka(i) .gt. 0.d0 .and. gammam(i) .gt. 0.d0 ) &
  dm(i) = soildepth(i) * gammam(i)
enddo

! if ksg(i) = 0.d0 -> no gw calculation
gw_switch = 0
do i = 1, num_of_landuse
 if( ksg(i) .gt. 0.d0 ) then
  gw_switch = 1
 else
  gammag(i) = 0.d0
  kg0(i) = 0.d0
  fpg(i) = 0.d0
  rgl(i) = 0.d0
 endif
enddo

end subroutine RRI_Read
