!*********************************************************************
subroutine Tecout_alloc(nx,ny,ValNum_temp)
    use tecout_mod
    implicit none
!---define
    integer :: nx,ny,ValNum_temp
!---
    call alloc_Vals(nx,ny,ValNum_temp)
    open(6, file = tecfile, status='unknown')
end subroutine Tecout_alloc
!*********************************************************************
subroutine Tecout_mkGrid(DX,DY,Z_temp)
    use tecout_mod
    implicit none
!---define
    integer :: i,j
    real(8) :: DX,DY,divN
    real(8) :: Y_temp(iMX),Z_temp(iMX-1,jMX-1)
!Xç¿ïWçÏê¨
      X(1)=0.
      do j=2,jMX
          X(j)=X(j-1)+DX
      end do
!Yç¿ïWçÏê¨
      Y_temp(1)=0.
      do i=2,iMX
          Y_temp(i)=Y_temp(i-1)+DY
      end do
      do i=iMX,1,-1
          Y(iMX-i+1)=Y_temp(i)
      end do
!Zç¿ïWçÏê¨-Cnt2Node
      do i=1,iMX-1
          do j=1,jMX-1
              Z_buf(i,j) = Z_temp(i,j)
              if (Z_buf(i,j)<0.) Z_buf(i,j)=0.
          end do
      end do
      do i=1,iMX
          do j=1,jMX
              if (i==1) then
                  if (j==1) then
                      if (Z_buf(i,j)>0.) then
                          Z(i,j)=Z_buf(i,j)
                      else
                          Z(i,j)=0.
                      end if
                  elseif (j==jMX) then
                      if (Z_buf(i,j-1)>0.) then
                          Z(i,j)=Z_buf(i,j-1)
                      else
                          Z(i,j)=0.
                      end if
                  else
                      divN=0.
                      if (Z_buf(i,j-1)>0.) divN=divN+1.
                      if (Z_buf(i,j  )>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)=(Z_buf(i,j-1)+Z_buf(i,j))/divN
                      else
                          Z(i,j)=0.
                      end if
                  end if
              elseif (i==iMX) then
                  if (j==1) then
                      if (Z_buf(i-1,j)>0.) then
                          Z(i,j)=Z_buf(i-1,j)
                      else
                          Z(i,j)=0.
                      end if
                  elseif (j==jMX) then
                      if (Z_buf(i-1,j-1)>0.) then
                          Z(i,j)=Z_buf(i-1,j-1)
                      else
                          Z(i,j)=0.
                      end if
                  else
                      divN=0.
                      if (Z_buf(i-1,j-1)>0.) divN=divN+1.
                      if (Z_buf(i-1,j  )>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)=(Z_buf(i-1,j-1)+Z_buf(i-1,j))/divN
                      else
                          Z(i,j)=0.
                      end if
                  end if
              else
                  if (j==1) then
                      divN=0.
                      if (Z_buf(i-1,j)>0.) divN=divN+1.
                      if (Z_buf(i  ,j)>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)=(Z_buf(i-1,j)+Z_buf(i,j))/divN
                      else
                          Z(i,j)=0.
                      end if
                  elseif (j==jMX) then
                      divN=0.
                      if (Z_buf(i-1,j-1)>0.) divN=divN+1.
                      if (Z_buf(i  ,j-1)>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)=(Z_buf(i-1,j-1)+Z_buf(i,j-1))/divN
                      else
                          Z(i,j)=0.
                      end if
                  else
                      divN=0.
                      if (Z_buf(i-1,j-1)>0.) divN=divN+1.
                      if (Z_buf(i  ,j-1)>0.) divN=divN+1.
                      if (Z_buf(i-1,j  )>0.) divN=divN+1.
                      if (Z_buf(i  ,j  )>0.) divN=divN+1.
                      if (divN>0.) then
                          Z(i,j)= &
                           ( Z_buf(i-1,j-1) + Z_buf(i  ,j  ) + &
                             Z_buf(i-1,j  ) + Z_buf(i  ,j-1) ) / divN
                      else
                          Z(i,j)=0.
                      end if
                  end if
              end if
          end do
      end do
end subroutine  Tecout_mkGrid
!*********************************************************************
subroutine Tecout_write_initialize( k, val1, val2, val3, val4)
    use Tecout_Mod
    implicit none
!---define
    integer :: i,j,k
    real(8) :: val1(iMX-1,jMX-1), val2(iMX-1,jMX-1), val3(iMX-1,jMX-1), val4(iMX-1,jMX-1)

    !ç≈èâÇ…ã§í ïîï™ÇèoóÕ
        write(6,'(a)') &
          ' VARIABLES = "X","Y","Z","width","depth","height","area_ratio" &
                        "Rain","River_H","River_Q","Surface_H","Surface_H_Max" '
        write(6,'(a,i5,a)') 'ZONE T="',k,'"'
        write(6,'(2(a,i5),a)') &
          'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK'
        write(6,'(a,i10)') 'SOLUTIONTIME=',k
        write(6,'(a)') 'VARLOCATION=([4-12]=CELLCENTERED)'
        
         do i=1,iMX
             do j=1,jMX
                 write(6,'(f12.3,",",$)') X(j)
             end do
             write(6,*)
         end do
         do i=1,iMX
             do j=1,jMX
                 write(6,'(f12.3,",",$)') Y(i)
             end do
             write(6,*)
         end do
         do i=1,iMX
             do j=1,jMX
                 write(6,'(f12.3,",",$)') Z(i,j)
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val1(i,j)
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val2(i,j)
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val3(i,j)
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val4(i,j)
             end do
             write(6,*)
         end do    
end subroutine Tecout_write_initialize    
!*********************************************************************
subroutine Tecout_write( k, val1, val2, val3, val4, area )
    use Tecout_Mod
    implicit none
!---define
    integer :: i,j,k
    real(8) :: val1(iMX-1,jMX-1), val2(iMX-1,jMX-1), val3(iMX-1,jMX-1), val4(iMX-1,jMX-1), area
                !"Rain","River_H","River_Q","Surface_H"
    if (k==1) then
    !ÇPâÒñ⁄
!        write(6,'(a)') &
!          'VARIABLES = "X","Y","Z","Rain","River_H","River_Q",Surface'
!        write(6,'(a,i5,a)') 'ZONE T="',k,'"'
!        write(6,'(2(a,i5),a)') &
!          'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK'
!        write(6,'(a,i10)') 'SOLUTIONTIME=',k
!        write(6,'(a)') 'VARLOCATION=([4-7]=CELLCENTERED)'
!        
!         do i=1,iMX
!         !do i=iMX,1,-1
!             do j=1,jMX
!                 write(6,'(f12.3,",",$)') X(j)
!             end do
!             write(6,*)
!         end do
!         do i=1,iMX
!         !do i=iMX,1,-1
!             do j=1,jMX
!                 write(6,'(f12.3,",",$)') Y(i)
!             end do
!             write(6,*)
!         end do
!         do i=1,iMX
!         !do i=iMX,1,-1
!             do j=1,jMX
!                 write(6,'(f12.3,",",$)') Z(i,j)
!             end do
!             write(6,*)
!         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val1(i,j) * 3600.d0 * 1000.d0
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val2(i,j)
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val3(i,j) * area
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val4(i,j)
             end do
             write(6,*)
         end do
         where(SufHmax < val4 ) SufHmax = val4
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') SufHmax(i,j)
             end do
             write(6,*)
         end do
         
      else
      !2âÒñ⁄à»ç~
        !write(6,'(a)') &
        !  'VARIABLES = "X","Y","Z","Rain","River_H","River_Q",Surface'
        write(6,'(a,i5,a)') 'ZONE T="',k,'"'
        write(6,'(2(a,i5),a)') &
          'I=' ,jMX, ',J=' ,iMX, ',DATAPACKING=BLOCK'
        write(6,'(a,i10)') 'SOLUTIONTIME=',k
        write(6,'(a)') 'VARLOCATION=([4-12]=CELLCENTERED)'
        write(6,'(a)') 'VARSHARELIST=([1-7]=1)'
        
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val1(i,j) * 3600.d0 * 1000.d0
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val2(i,j)
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val3(i,j) * area
             end do
             write(6,*)
         end do
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') val4(i,j)
             end do
             write(6,*)
         end do
         where(SufHmax < val4 ) SufHmax = val4
         do i=1,iMX-1
             do j=1,jMX-1
                 write(6,'(f12.5,",",$)') SufHmax(i,j)
             end do
             write(6,*)
         end do
         
      end if
end subroutine Tecout_write


