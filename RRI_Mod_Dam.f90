module dam_mod

integer dam_switch
character*256 damfile                            ! dam control file
    
integer :: dam_num                               ! number of dam
character*256, allocatable :: dam_name(:)        ! dam name

integer, allocatable :: dam_kind(:)              ! dam id number
integer, allocatable :: dam_ix(:), dam_iy(:)     ! dam location (x, y)
integer, allocatable :: dam_loc(:)               ! dam location (k)
integer, allocatable :: dam_state(:)             ! dam status (full:1, not full:0)
integer, allocatable :: damflg(:)                ! dam exist at each grid-cell (exist:1, not exist:0)

real(8), allocatable :: dam_qin(:)               ! dam inflow
real(8), allocatable :: dam_vol(:)               ! storage volume
real(8), allocatable :: dam_vol_temp(:)          ! storage volume (temporary)
real(8), allocatable :: dam_volmax(:)            ! maximum storage
real(8), allocatable :: dam_qout(:)              ! dam outflow
real(8), allocatable :: dam_floodq(:)            ! flood discharge
real(8), allocatable :: dam_maxfloodq(:)         ! max flood discharge (for non-const. cont.)
real(8), allocatable :: dam_rate(:)              ! discharge increasing rate

end module dam_mod
