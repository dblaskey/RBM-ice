SUBROUTINE Water_Balance
!
USE Block_Hydro
USE Block_Network
!
IMPLICIT NONE
!
integer     :: nc,ncnf,nr,nrc,ntrib,nntrib,total_cells
integer     :: tot_trib_in,trib_seg_in
real        :: Q_start,Q_sum_trib
!
total_cells=0
!
do nr = 1,nreach
! 
  do nc = 1,no_cells(nr)
    total_cells=total_cells+1
    nrc = reach_cell(nr,nc)
    if (nc .eq. 1) then
!
!  Flow into the headwaters cell is the same as the flow out
!
      Q_in(nrc) = Q_out(nrc)
    else
!
!  Flow into downstream cells is equal to the flow out of the most upstream cell
!
      Q_in(nrc) = Q_out(nrc-1)
    end if
    ntrib=no_tribs(nrc)
    Q_sum_trib = 0.0
    if (ntrib .gt. 0) then
      do nntrib  = 1,ntrib
        ncnf = conflnce(nrc,nntrib)
        Q_sum_trib = Q_sum_trib + Q_trib(ncnf)
      end do
    end if
!
! Distributed flow is divided equally among each the segments in a cell
!
    Q_diff(nrc) = (Q_out(nrc) - Q_in(nrc) - Q_sum_trib)/ndelta(nrc)
!
!    
  end do
end do
end SUBROUTINE Water_Balance