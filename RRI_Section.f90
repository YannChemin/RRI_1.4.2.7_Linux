! RRI_Section

subroutine set_section
use globals
implicit none

integer i, j, k, id, div_max
real(8) rdummy
character*256 sec_file_name
character*6 id_char

allocate( sec_div(sec_id_max), sec_depth(sec_id_max), sec_height(sec_id_max) )

div_max = 0
do i = 1, ny
 do j = 1, nx

  id = sec_map(i, j)
  if( id .gt. 0 ) then
   call int2char( id, id_char )
   sec_file_name = trim(sec_file) // trim(id_char) // ".txt"

   open(10, file = sec_file_name, status = 'old')
   read(10, *) sec_div(id), sec_depth(id), sec_height(id)
   
   if(div_max .le. sec_div(id)) div_max = sec_div(id)
   close(10)
  endif

 enddo
enddo

allocate( sec_hr(sec_id_max, div_max), sec_area(sec_id_max, div_max) )
allocate( sec_peri(sec_id_max, div_max), sec_b(sec_id_max, div_max) )
allocate( sec_ns_river(sec_id_max, div_max) )

do i = 1, ny
 do j = 1, nx

  id = sec_map(i, j)
  if( id .gt. 0 ) then
   call int2char( id, id_char )
   sec_file_name = trim(sec_file) // trim(id_char) // ".txt"

   open(10, file = sec_file_name, status = 'old')
   read(10, *) sec_div(id), sec_depth(id), sec_height(id)

   depth(i, j) = sec_depth(id)
   height(i, j) = sec_height(id)
   if(height(i, j) .lt. 0.d0) height(i, j) = 0.d0 ! added ver 1.4.2.4

   do k = 1, sec_div(id)
    read(10, *) sec_hr(id, k), sec_peri(id, k), sec_b(id, k), sec_ns_river(id, k)
   enddo

   sec_area(id, 1) = sec_b(id, 1) * sec_hr(id, 1)
   do k = 2, sec_div(id)
    sec_area(id, k) = sec_area(id, k-1) + sec_b(id, k) * (sec_hr(id, k) - sec_hr(id, k-1))
   enddo

   !width(i, j) = sec_b(id, div_max)
   width(i, j) = sec_b(id, sec_div(id)) ! modified on Feb 14, 2021

   close(10)

   riv(i, j) = 1 ! added by T.Sayama on Dec. 27, 2021
   
  endif
 enddo
enddo

return
end


subroutine sec_hq_riv(h, dh, k, q)
use globals
implicit none

real(8) h, dh, q, p, n, a, b
integer k, id, div_max, i

id = sec_map_idx(k)
div_max = sec_div(id)

if( h .le. sec_hr(id, 1) ) then

 p = sec_peri(id, 1)
 n = sec_ns_river(id, 1)
 a = h * sec_b(id, 1)

elseif( h .gt. sec_hr(id, div_max) ) then

 p = sec_peri(id, div_max)
 n = sec_ns_river(id, div_max)
 a = sec_area(id, div_max) + (h - sec_hr(id, div_max)) * sec_b(id, div_max)

else

 do i = 2, div_max

  if( h .le. sec_hr(id, i) ) then
   p = sec_peri(id, i)
   n = sec_ns_river(id, i)
   a = sec_area(id, i - 1) + (h - sec_hr(id, i - 1)) * sec_b(id, i)
   exit
  endif

 enddo
endif

q = 1.d0 / n  * ( a / p ) ** (2.d0 / 3.d0) * sqrt( abs(dh) ) * a

return
end


subroutine hr2vr(hr, k, vr)
use globals
implicit none

real(8) hr, vr, a
integer k, id, div_max, i

id = sec_map_idx(k)

if( id .le. 0 ) then

 vr = hr * area * area_ratio_idx(k)

else

 div_max = sec_div(id)

 if( hr .le. sec_hr(id, 1) ) then

  a = hr * sec_b(id, 1)

 elseif( hr .gt. sec_hr(id, div_max) ) then

  a = sec_area(id, div_max) + (hr - sec_hr(id, div_max)) * sec_b(id, div_max)

 else

  do i = 2, div_max

   if( hr .le. sec_hr(id, i) ) then
    a = sec_area(id, i - 1) + (hr - sec_hr(id, i - 1)) * sec_b(id, i)
    exit
   endif

  enddo
 endif

 vr = a * len_riv_idx(k)

endif

return
end


subroutine vr2hr( vr, k, hr )
use globals
implicit none

real(8) hr, vr, a
integer k, id, div_max, i

id = sec_map_idx(k)

if( id .le. 0 ) then

 hr = vr / ( area * area_ratio_idx(k) )

else

 div_max = sec_div(id)

 a = vr / len_riv_idx(k)

 if( a .le. sec_area(id, 1) ) then

  hr = a / sec_b(id, 1)

 elseif( a .gt. sec_area(id, div_max) ) then

  hr = (a - sec_area(id, div_max)) / sec_b(id, div_max) + sec_hr(id, div_max)

 else

  do i = 2, div_max

   if( a .le. sec_area(id, i) ) then
    hr = (a - sec_area(id, i-1)) / sec_b(id, i) + sec_hr(id, i - 1)
    exit
   endif

  enddo
 endif

endif

return
end


subroutine hr_update(hr_org, vr_inc, k, hr_new)
use globals
implicit none

real(8) hr_org, vr_inc, hr_new, vr_org, vr_new
integer k

call hr2vr(hr_org, k, vr_org)
vr_new = vr_org + vr_inc
call vr2hr(vr_new, k, hr_new)

return
end


subroutine sec_h2b(h, k, b)
use globals
implicit none

real(8) h, b
integer k, id, div_max, i

id = sec_map_idx(k)
if( id .le. 0 ) then

 b = width_idx(k)

else

 div_max = sec_div(id)

 do i = 1, div_max
  if( h .le. sec_hr(id, i) ) then
   b = sec_b(id, i)
   exit
  endif
  b = sec_b(id, div_max)
 enddo

endif

return
end
