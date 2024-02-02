! RRI_Break.f90

! current version of RRI model does not call this subroutine
! this subroutine is still under development

!****** Levee Break Model ***************************

subroutine levee_break(t, hr, hs)
use globals

implicit none

real(8) a, b, vn, p, hr(ny, nx), hs(ny, nx), hw
real(8) x, y
real ran
integer t, every_t, i, j

integer seed_param, se
parameter( seed_param = 51234592 )

every_t = 144 ! -> every day (when dt = 600 sec)

if( mod(t, every_t) .ne. 0 ) return

! first time
if( t .eq. every_t ) then

 !se = seed_param
 call system_clock(se)

 call seed(se)
 open(301, file = "levee_break.out")
 write(301,'("id, x, y, t, hr, ran, p")')
 id_break = 0
endif

! logit model parameter
! hc = 10m
! hw = 7m -> p = 0.2
! hw = 10m -> p = 0.99
a = 18.55
b = -13.96

do i = 1, ny
 do j = 1, nx

  if( domain(i,j) .eq. 0 ) cycle
  !if( riv(i,j) .eq. 0 ) cycle
  if( height(i, j) .le. 0 ) cycle

  if( riv(i, j) .eq. 1 ) then
   ! levee on river grid cells
   hw = hr(i, j) - depth(i, j)
   if( hw .le. 0.1d0 ) cycle
  else
   ! levee on slope grid cells
   hw = 0.d0
   if( j.ne.nx ) hw = max( hw, hs(i,j+1) )
   if( i.ne.ny .or. j.ne.nx ) hw = max( hw, hs(i+1,j+1) )
   if( i.ne.ny ) hw = max( hw, hs(i+1,j) )
   if( i.ne.ny .or. j.ne.1 ) hw = max( hw, hs(i+1,j-1) )
   if( j.ne.1 ) hw = max( hw, hs(i,j-1) )
   if( i.eq.1 .or. j.eq.1 )hw = max( hw, hs(i-1,j-1) )
   if( i.eq.1 ) hw = max( hw, hs(i-1,j) )
   if( i.eq.1 .or. j.eq.1) hw = max( hw, hs(i-1,j-1) )
  endif

  vn = a + b * height(i, j) / hw
  p = 1 / (1 + exp( - vn )) ! break probability for hw

  call random(ran)
  !if( ran .lt. 0.1 ) ran = 0.1

  if( ran .lt. p ) then
   ! break
   id_break = id_break + 1
   x = xllcorner + (j - 0.5) * cellsize
   y = yllcorner + (ny - i + 0.5) * cellsize

   ! levee on slope grid cells
   if( riv(i, j) .ne. 1 ) then
    zb_slo_idx(slo_ij2idx(i, j)) = zb_slo_idx(slo_ij2idx(i, j)) - height(i, j)
   endif

   height(i, j) = 0.d0

   write(*,*) "Levee Break ", id_break, acc(i, j), hw
   write(301,'(i4, ",", f12.3, ",", f12.3, ",", i7, ",", f12.3, ",", f12.3, ",", f12.3)') id_break, x, y, t, hw, ran, p
  endif

 enddo
enddo

return
end
