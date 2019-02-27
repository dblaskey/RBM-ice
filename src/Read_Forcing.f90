SUBROUTINE Read_Forcing
!
USE Block_Energy
USE Block_Hydro
USE Block_Network
! 
IMPLICIT NONE
!
integer:: nc,ncell,nnd,no_flow,no_heat,nr,nrec_flow,nrec_heat
real:: Q_avg


no_flow=0
no_heat=0
do nr=1,nreach
  do nc=1,no_cells(nr)-1 
    no_flow=no_flow+1
    no_heat=no_heat+1
!
    nrec_flow=flow_cells*(ndays-1)+no_flow
    nrec_heat=heat_cells*(ndays-1)+no_heat
!
    read(35,'(2i5,2f10.1,2f6.1,f7.1,f6.2)' &
           ,rec=nrec_flow) nnd,ncell &
           ,Q_in(no_heat),Q_out(no_heat),Q_diff(no_heat) &  
           ,depth(no_heat),width(no_heat),u(no_heat)
!
! Constant flow scenario
!
    u(no_heat) = 2.5
!    Q_in(no_heat)  = 1000.
!    Q_out(no_heat) = 1000.
!    Q_out(15) = 2000.
!    if (no_heat .ge. 16) then
!      Q_in(no_heat)  = 2000.
!      Q_out(no_heat) = 2000.
!    end if
!    if(u(no_heat).lt.0.5) u(no_heat)=0.5
!
    if(ncell.ne.no_heat) write(*,*) 'Flow file error',ncell,no_heat 
!
    read(36,'(i5,2f6.1,2f7.4,f6.3,f7.1,f5.1)' &
           ,rec=nrec_heat) ncell &
           ,dbt(no_heat),ea(no_heat) &
           ,Q_ns(no_heat),Q_na(no_heat),rho &
           ,press(no_heat),wind(no_heat)
!          
!  if(ncell.ne.no_heat) write(*,*) 'Heat file error',ncell,no_heat 
!
!  Added variable ndelta (UW_JRY_2011/03/15

!

    delta_n=ndelta(ncell)
! 
    Q_avg=0.5*(Q_in(no_heat)+Q_out(no_heat))
!
    dt(no_heat)=dx(no_heat)/u(no_heat)

!
!  Added check to see if travel time of parcel exceeds the
!  computational interval.  If so, it writes to file fort.45.
!  One should check to see if there are NaN's at that node.
!  (UW_JRY_2011/03/15)
!
    if(dt(no_heat).gt.dt_comp) write(45,*) &
           'Travel time=',dt(no_heat) &
            , '> dt_comp at node -',no_heat
  end do

!
!       Read the meteorology for the last cell, but not the flow
!
  no_heat=no_heat+1 
  Q_in(no_heat)=Q_out(no_heat-1)
  Q_out(no_heat)=Q_in(no_heat)
  u(no_heat)=u(no_heat-1)
  depth(no_heat)=depth(no_heat-1)
  width(no_heat)=width(no_heat-1)
  dt(no_heat)=dx(ncell)/u(no_heat)
!
!      Tributary flow is Q_out from the next to the last cell
!
  Q_trib(nr)=Q_out(no_heat)    
!
!  Meteorology for the part of the cell that's in the confluent
!
  nrec_heat=heat_cells*(ndays-1)+no_heat
  read(36,'(i5,2f6.1,2f7.4,f6.3,f7.1,f5.1)' &
         ,rec=nrec_heat) ncell &
         ,dbt(no_heat),ea(no_heat) &   
         ,Q_ns(no_heat),Q_na(no_heat),rho &
         ,press(no_heat),wind(no_heat)
end do
!  
!  Call the Water Balance routine
!
  call Water_Balance
!
END SUBROUTINE Read_Forcing
