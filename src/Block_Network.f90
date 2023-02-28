Module Block_Network
!
! Module with stream topology variables
!
    integer, dimension(:), allocatable  ::no_celm,no_cells,no_tribs
    integer, dimension(:), allocatable  ::head_cell
!
    integer, dimension(:,:), allocatable::conflnce,reach_cell,segment_cell,trib
!
!
! Integer variables 
!
    integer             :: flow_cells,heat_cells
    integer             :: n1,n2,ndays,nreach,ntrb,nwpd
    integer,parameter   :: ns_max=3000
    integer             :: start_year,start_month,start_day
    integer             :: end_year,end_month,end_day
!
! Real variables
!
    real, parameter                   :: n_default=2
    real                              :: delta_n
    real                              :: dt_comp,x_bndry
    real, dimension(:), allocatable   :: ndelta
!
end module Block_Network
