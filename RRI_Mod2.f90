! Variables used for "Adaptive Stepsize Control Runge-Kutta"

module runge_mod

real(8) eps, ddt_min_riv, ddt_min_slo
! for standard simulation
parameter( eps = 0.01d0 )
parameter( ddt_min_riv = 0.1d0 )
parameter( ddt_min_slo = 1.d0 )

! for detailed simulation (new setting for detail)
!parameter( eps = 0.005d0 )
!parameter( ddt_min_riv = 0.1d0 )
!parameter( ddt_min_slo = 1.d0 )

! for rough simulation
!parameter( eps = 0.1d0 )
!parameter( ddt_min_riv = 1.d0 )
!parameter( ddt_min_slo = 1.d0 )

real(8) safety, pgrow, pshrnk, errcon
parameter (safety=0.9d0,pgrow=-.2d0,pshrnk=-.25d0,errcon=1.89d-4)

real(8), allocatable, save :: ks2(:), ks3(:), ks4(:), ks5(:), ks6(:)
real(8), allocatable, save :: hs_temp(:), hs_err(:)

real(8), allocatable, save :: kg2(:), kg3(:), kg4(:), kg5(:), kg6(:)
real(8), allocatable, save :: hg_temp(:), hg_err(:)

real(8), allocatable, save :: kr2(:), kr3(:), kr4(:), kr5(:), kr6(:)
!real(8), allocatable, save :: hr_temp(:), hr_err(:)
real(8), allocatable, save :: vr_temp(:), vr_err(:), hr_err(:) ! v1.4

real(8) a2,a3,a4,a5,a6,b21,b31,b32,b41,b42,b43,b51,b52,b53, &
        b54,b61,b62,b63,b64,b65,c1,c3,c4,c6,dc1,dc3,dc4,dc5,dc6
parameter (a2=.2d0,a3=.3d0,a4=.6d0,a5=1.d0,a6=.875d0,b21=.2d0,b31=3.d0/40.d0, &
     b32=9.d0/40.d0,b41=.3d0,b42=-.9d0,b43=1.2d0,b51=-11.d0/54.d0,b52=2.5d0, &
     b53=-70.d0/27.d0,b54=35.d0/27.d0,b61=1631.d0/55296.d0,b62=175.d0/512.d0, &
     b63=575.d0/13824.d0,b64=44275.d0/110592.d0,b65=253.d0/4096.d0,c1=37.d0/378.d0, &
     c3=250.d0/621.d0,c4=125.d0/594.d0,c6=512.d0/1771.d0,dc1=c1-2825.d0/27648.d0, &
     dc3=c3-18575.d0/48384.d0,dc4=c4-13525.d0/55296.d0,dc5=-277.d0/14336.d0, &
     dc6=c6-.25d0)

real(8), save :: errmax

end module runge_mod
