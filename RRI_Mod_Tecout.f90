module tecout_mod
    integer :: tec_switch
    character*256 tecfile
    integer :: iMX,jMX,ValNum
    real(8),allocatable :: X(:),Y(:),Z(:,:),Z_buf(:,:)
    real(8),allocatable :: SufHmax(:,:)
    !--------
    contains
    !--------
    subroutine alloc_Vals(nx,ny,ValNum_temp)
        integer :: nx,ny,ValNum_temp
        iMX = ny + 1
        jMX = nx + 1
        ValNum = ValNum_temp
        allocate( X(jMX),Y(iMX),Z(iMX,jMX),Z_buf(ny,nx) )
        allocate( SufHmax(iMX-1,jMX-1) )
        SufHmax(:,:) = -0.1d0
    end subroutine alloc_Vals
end module tecout_mod
