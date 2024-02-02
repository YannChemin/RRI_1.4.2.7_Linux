! RRI_Dam.f90
! modified by T.Sayama on Dec. 2, 2020

! reading dam control file
subroutine dam_read
use globals
use dam_mod
implicit none
    
integer :: i, ios
    
allocate( damflg(riv_count), dam_qin(riv_count) )
damflg(:) = 0

if( dam_switch.eq.1 ) then

 open(99, file=damfile, status="old" )
 read(99,*) dam_num
 allocate( dam_name(dam_num), dam_kind(dam_num) &
           , dam_ix(dam_num), dam_iy(dam_num) &
           , dam_vol(dam_num), dam_vol_temp(dam_num) &
           , dam_volmax(dam_num), dam_state(dam_num) &
           , dam_qout(dam_num), dam_loc(dam_num) &
           , dam_floodq(dam_num), dam_maxfloodq(dam_num) &
           , dam_rate(dam_num))

 dam_vol(:) = 0.d0
 dam_state(:) = 0
 dam_floodq(:) = 0.d0
 dam_maxfloodq(:) = 0.d0
 dam_rate(:) = 0.d0

 read(99,*,iostat = ios) dam_name(1), dam_iy(1), dam_ix(1), dam_volmax(1), dam_floodq(1), dam_maxfloodq(1), dam_rate(1), dam_vol(1)
 rewind(99)
 read(99,*) dam_num
 !if(ios .gt. 0) then ! old version (not recommended) bug fix on Nov 27, 2021
 if(ios .ne. 0) then ! old version (not recommended)
  do i = 1, dam_num
   read(99,*) dam_name(i), dam_iy(i), dam_ix(i), dam_volmax(i), dam_floodq(i)
   dam_loc(i) = riv_ij2idx( dam_iy(i), dam_ix(i) )
   damflg(dam_loc(i)) = i
  enddo
 else ! new version
  do i = 1, dam_num
   read(99,*) dam_name(i), dam_iy(i), dam_ix(i), dam_volmax(i), dam_floodq(i), dam_maxfloodq(i), dam_rate(i), dam_vol(i)
   dam_loc(i) = riv_ij2idx( dam_iy(i), dam_ix(i) )
   damflg(dam_loc(i)) = i
  enddo
 end if
 close(99)
end if
end subroutine dam_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! calculating inflow to dam
subroutine dam_prepare(qr_idx)

use globals
!use dam_mod, only :dam_qin
use dam_mod
implicit none

integer :: i, k, kk
real(8) :: qr_idx(riv_count), vr_idx(riv_count)

!dam_qin(:) = 0.d0
!do k = 1, riv_count
! kk = down_riv_idx(k)
! dam_qin(kk) = dam_qin(kk) + qr_idx(k)
!enddo

! modified by TS on June 16, 2016
dam_qin(:) = qr_idx(:)

end subroutine dam_prepare

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dam_operation(k)

use globals, only :ddt
use dam_mod
implicit none

integer :: k
real(8) :: qdiff
    
dam_qout(damflg(k)) = 0.d0

if ( dam_maxfloodq(damflg(k)) .eq. 0 .and. dam_rate(damflg(k)) .eq. 0 ) then
 ! constant flood peak cut
 if ( dam_qin(k) .lt. dam_floodq(damflg(k)) ) then
  if( dam_vol(damflg(k)) .le. 0 ) then ! ’Ç‰Á
   dam_qout(damflg(k)) = dam_qin(k)
  else ! ’Ç‰Á
   if( dam_qin(k) .lt. 0.25 * dam_floodq(damflg(k)) ) then
    dam_qout(damflg(k)) = 0.25 * dam_floodq(damflg(k)) ! ’Ç‰Á
    qdiff = (dam_qin(k) - dam_qout(damflg(k))) ! ’Ç‰Á
    dam_vol_temp(damflg(k)) = dam_vol_temp(damflg(k)) + qdiff * ddt ! ’Ç‰Á
   else
    dam_qout(damflg(k)) = dam_qin(k)
   endif
  endif
 else
  if (dam_state(damflg(k)) .eq. 0) then
  ! still have space
   dam_qout(damflg(k)) = dam_floodq(damflg(k))
   qdiff = (dam_qin(k) - dam_qout(damflg(k)))
   dam_vol_temp(damflg(k)) = dam_vol_temp(damflg(k)) + qdiff * ddt
  else
   ! no more space
   dam_qout(damflg(k)) = dam_qin(k)
  end if
 end if

else
 ! non-constant flood peak cut
 if ( dam_qin(k) .lt. dam_floodq(damflg(k)) ) then
  if( dam_vol(damflg(k)) .le. 0 ) then ! ’Ç‰Á
   dam_qout(damflg(k)) = dam_qin(k)
  else ! ’Ç‰Á
   if( dam_qin(k) .lt. 0.25 * dam_floodq(damflg(k)) ) then
    dam_qout(damflg(k)) = 0.25 * dam_floodq(damflg(k)) ! ’Ç‰Á
    qdiff = (dam_qin(k) - dam_qout(damflg(k))) ! ’Ç‰Á
    dam_vol_temp(damflg(k)) = dam_vol_temp(damflg(k)) + qdiff * ddt ! ’Ç‰Á
   else
    dam_qout(damflg(k)) = dam_qin(k)
   endif
  endif
 else
  if (dam_state(damflg(k)) .eq. 0) then
  ! still have space
   dam_qout(damflg(k)) = dam_floodq(damflg(k)) + &
      dam_rate(damflg(k)) * (dam_qin(k) - dam_floodq(damflg(k)))
   if( dam_qout(damflg(k)) .gt. dam_maxfloodq(damflg(k)) ) &
      dam_qout(damflg(k))  = dam_maxfloodq(damflg(k))
   qdiff = (dam_qin(k) - dam_qout(damflg(k)))
   dam_vol_temp(damflg(k)) = dam_vol_temp(damflg(k)) + qdiff * ddt
  else
   ! no more space
   dam_qout(damflg(k)) = dam_qin(k)
  end if
 end if
end if

end subroutine dam_operation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gate_operation(k, hr_k, hr_kk)

use globals, only :down_riv_idx, zb_riv_idx, depth, riv_idx2i, riv_idx2j
use dam_mod
implicit none

integer :: k, kk, i, j
real(8) :: hr_k, hr_kk

dam_qout(damflg(k)) = 0.d0
kk = down_riv_idx(k)
!i = riv_idx2i(k)
!j = riv_idx2j(k)
dam_qout(damflg(k)) = dam_qin(k)
if( dam_qout(damflg(k)) .ge. dam_floodq(damflg(k)) ) &
 dam_qout(damflg(k)) = dam_floodq(damflg(k))
! zb_riv = zs - depth
! dam_floodq -> gate_close_level
!if ( hr_kk + zb_riv_idx(kk) - zb_riv_idx(k) .gt. dam_floodq(damflg(k)) ) then
! dam_qout(damflg(k)) = 0.d0 ! gate close
 ! dam_floodq -> pump_operate_level
!if( (hr_k - depth(i, j)) .gt. dam_maxfloodq(damflg(k)) ) then
! dam_qout(damflg(k)) = dam_rate(damflg(k))
! if( dam_qout(damflg(k)) .ge. dam_qin(k) ) dam_qout(damflg(k)) = dam_qin(k)
!endif
!else
! dam_qout(damflg(k)) = dam_qin(k)
!endif

end subroutine gate_operation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dam_checkstate(qr_ave_idx)

use globals,only :riv_count,area
use dam_mod
implicit none
    
real(8) :: qr_ave_idx(riv_count)
integer :: i
    
do i=1, dam_num
 dam_vol_temp(i) = dam_vol_temp(i) / 6.d0
 dam_vol(i) = dam_vol(i) + dam_vol_temp(i)
 dam_state(i) = 0 ! added on January 8, 2021
 if (dam_vol(i) .gt. dam_volmax(i)) dam_state(i) = 1
end do
    
end subroutine dam_checkstate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dam_write

!use globals, only :area
use dam_mod
implicit none
integer i

write(1001,'(10000f20.5)') (dam_vol(i), dam_qin(dam_loc(i)), dam_qout(i), i = 1, dam_num) ! v1.4

end subroutine dam_write

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine dam_write_cnt

!use globals, only :area
use dam_mod
implicit none
integer i

write(1002,'(i10)') dam_num
do i = 1, dam_num
 write(1002,'(a10, 2i8, 10000f20.5)') dam_name(i), dam_iy(i), dam_ix(i), dam_volmax(i), dam_floodq(i), &
         dam_maxfloodq(i), dam_rate(i), dam_vol(i)
enddo

end subroutine dam_write_cnt
