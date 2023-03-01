SUBROUTINE SYSTMM(flowPrefix,heatPrefix,outPrefix)
!
use Block_Energy
use Block_Hydro
use Block_Ice_Snow
use Block_Network
!
Implicit None
! 
character (len=200):: temp_file
character (len=200):: outPrefix
character (len=200):: flowPrefix
character (len=200):: heatPrefix
character (len=200):: flow_file
character (len=200):: heat_file
character (len=200):: year_file_id
!
integer          :: ncell, nc_head
integer          :: nd, ndd, nr, ns
integer          :: nn, nobs, nyear, nd_year, ntmp
!
real             :: hpd, xd, xdd, xwpd, xd_year, year
real             :: q_rslt,T_mohseni,xn_avg
real(8)          :: time
!
! Allocate the arrays
!
allocate (SNOW(heat_cells))
allocate (SUB_ZERO(heat_cells))
SNOW= .FALSE.
SUB_ZERO = .FALSE.
!
allocate (T_head(nreach))
allocate (T_trib(nreach))
allocate (depth(heat_cells))
allocate (Q_in(heat_cells))
allocate (Q_out(heat_cells))
allocate (Q_diff(heat_cells))
allocate (Q_trib(nreach))
allocate (width(heat_cells))
allocate (u(heat_cells))
allocate (dt(2*heat_cells))
allocate (dbt(heat_cells))
allocate (ea(heat_cells))
allocate (Qns(heat_cells))
allocate (Qna(heat_cells))
allocate (press(heat_cells))
allocate (wind(heat_cells))
!
allocate(ICE(heat_cells))
ICE = 100.
allocate(ice_temp(nreach,-2:ns_max,2))
ice_temp = 0.001
allocate(ice_thick(nreach,-2:ns_max,2))
ice_thick = -0.001!
!
! Initialize some arrays
!
dt_part=0.
x_part=0.
no_dt=0
nstrt_elm=0
temp=0.5
! Initialize headwaters temperatures
!
T_head=0.1 ! Updated for Alaska by DJB 3/1/23
!
!
! Initialize smoothed air temperatures for estimating headwaters temperatures
!
xn_avg = nn_avg
T_smth=0.0
tmp_arry=0.0
!
!
!
! Initialize dummy counters that facilitate updating simulated values
!
n1=1
n2=2
nobs=0
ndays=0
xwpd=nwpd
hpd=1./xwpd
!
!     Year loop starts
!
do nyear=start_year,end_year
  write(*,*) ' Simulation Year - ',nyear,start_year,end_year
  !
  ! Open yearly output file (Added by DJB 3/1/23)
  write(year_file_id, '(i0)') nyear
  temp_file=TRIM(outPrefix)//'_'//TRIM(ADJUSTL(year_file_id))//'.temp'
  open(20,file=TRIM(temp_file),status='unknown')
  !
  !     Open file with yearly hydrologic data (Added by DJB 3/1/23)
  !
  flow_file=TRIM(flowPrefix)//'_'//TRIM(ADJUSTL(year_file_id))
  open(unit=35,FILE=TRIM(flow_file), ACCESS='SEQUENTIAL',FORM='FORMATTED', STATUS='old')
  !
  !     Open file with yearly meteorologic data (Added by DJB 3/1/23)
  !
  heat_file=TRIM(heatPrefix)//'_'//TRIM(ADJUSTL(year_file_id))
  open(unit=36,FILE=TRIM(heat_file), ACCESS='SEQUENTIAL',FORM='FORMATTED', STATUS='old')
  !
  nd_year=365
  if (mod(nyear,4).eq.0) nd_year=366
!
!     Day loop starts
!
  DO nd=1,nd_year
    year=nyear
    xd=nd
    xd_year=nd_year
!     Start the numbers of days-to-date counter
!
    ndays=ndays+1
!
!    Daily period loop starts
!
      DO ndd=1,nwpd
      xdd = ndd
      time=year+(xd+(xdd-0.5)*hpd)/xd_year 
!
! Read the hydrologic and meteorologic forcings
!
        call READ_FORCING(nyear,nd)
!
!     Begin reach computations
!
!
!     Begin cycling through the reaches
!
      do nr=1,nreach
!
        nc_head=segment_cell(nr,1)
!
!     Determine smoothing parameters (UW_JRY_2011/06/21)
!
!        rminsmooth=1.0-smooth_param(nr)
!        T_smth(nr)=rminsmooth*T_smth(nr)+smooth_param(nr)*dbt(nc_head)
      do nn = 2,nn_avg
        tmp_arry(nn) = T_smth(nr,nn-1)
        tmp_arry(1) = dbt(nr)
        T_Mohseni = sum(tmp_arry)/xn_avg
      end do
      do nn = 1,nn_avg
        T_smth(nr,nn) = tmp_arry(nn)
      end do

!     
!     Variable Mohseni parameters (UW_JRY_2011/06/16)
! 
        T_head(nr)=mu(nr)+(alphaMu(nr) &
                  /(1.+exp(gmma(nr)*(beta(nr)-T_Mohseni)))) 
!
      temp(nr,0,n1)=T_head(nr)
      temp(nr,1,n1)=T_head(nr)
      x_bndry=x_dist(nr,0) - 1.0
!
! Begin cell computational loop
!
        do ns=1,no_celm(nr)  
!
          ncell = segment_cell(nr,ns)
!
          if (ice_temp(nr,ns,n2) .gt. -0.01) then
            ICE(ncell) = 102
            ice_thick(nr,ns,n2) = 0.000
!           
             call Ice_Free (nd,nr,ns,ncell,nc_head)
          else
!
             call ENERGY_ICE (q_rslt, nd, ncell, nr, ns)
!
          end if
!
!
!   Write all temperature output UW_JRY_11/08/2013
!   The temperature is output at the beginning of the 
!   reach.  It is, of course, possible to get output at
!   other points by some additional code that keys on the
!   value of ndelta (now a vector)(UW_JRY_11/08/2013)
!
!        if (ice_thick(nr,ns,n2) .lt. 0.009) ICE(ncell) = 100.
        do nseg_temp=1,nseg_out_num
            if (nseg_out(nr,ncell,nseg_temp).eq.ns) then
                call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),depth(ncell), &
                   Q_in(ncell),ice_thick(nr,ns,n2),ICE(ncell))
            !if (ncell.eq.1827) write(89,*)nyear,nd,nr,ncell,ns,T_0
            end if
        end do
!
!        call WRITE(time,nd,nr,ncell,ns,T_0,T_head(nr),dbt(ncell),depth(ncell), &
!                   Q_in(ncell),ice_thick(nr,ns,n2),ICE(ncell))
!
!     End of computational element loop
!
              end do
!     End of reach loop
!
            end do
            ntmp=n1
            n1=n2
            n2=ntmp
!
!     End of weather period loop (NDD=1,NWPD)
!
          end do
!
!     End of main loop (ND=1,365/366)
!
        end do
        CLOSE(20) ! Close yearly output file (Added DJB 3/1/23)
        CLOSE(35) ! Close yearly flow file (Added DJB 3/1/23)
        CLOSE(36) ! Close yearly heat file (Added DJB 3/1/23)
!
!     End of year loop
!
      end do
!
!
!     ******************************************************
!                        return to rmain
!     ******************************************************
!
end SUBROUTINE SYSTMM
