Module Block_Ice_Snow
!
! Module with ice and snow parameters and variables
!
! Integer variables 
!
!
! Logical variables 
!
    logical,dimension(:),allocatable :: h_rate,SNOW, SUB_ZERO
!
! Real variables
!
    real, parameter                 :: ice_albedo = 0.25
    real, parameter                 :: ice_albedo = 0.8
    real,parameter                  :: alpha_ice = 0.5      ! Ice reflectivity - Parkinson-Wash
    real,parameter                  :: Ch_Cg = 1.75e-03     ! Transfer coefficent - Parkinson-Wash
    real,parameter                  :: kappa_ice = 8.0      ! J/sec/m**2/DegK - Wanders et al
    real, parameter                 :: ice_cndctvy = 1.04   ! W/m/deg K-Parkinson 
    real, parameter                  :: h2o_albedo = 0.05 
    real, parameter                 :: h2o_dnsty = 1000.0   ! Kg/m**3
    real, parameter                 :: h_ice_min = 0.025    ! meters
    real, parameter                 :: I0 = 0.17            ! Fraction of SW thru ice surface
    real, parameter                 :: ice_min_thck = 0.001 ! meters
    real, parameter                 :: sno_dnsty = 330.0    ! Kg/m**3
    real, parameter                 :: sens_h2o_ice = 0.280e-03 !Saucier et al
    real, parameter                 :: T_Kelvin = 273.15    ! Freezing point of H2O, degK
    real, parameter                 :: ltnt_ice_air = 0.840e-03 !Saucier et al
!
    real, dimension(:),allocatable  :: ICE
!
    real, dimension(:,:,:), allocatable :: ice_temp, ice_thick, snow_temp
!
!    real, parameter, dimension(2)   :: snow_albedo = (/0.85,0.60/)
    real, parameter, dimension(2)   :: ice_extnct  = (/1.5,20.0/) ! 1/meters
    real, parameter, dimension(2)   :: snow_extnct = (/6.0,20.0/) ! 1/meters
    real,dimension(2),parameter     :: a_vapor = (/9.5,7.5/)    !Parkinson and Washington
    real,dimension(2),parameter     :: b_vapor = (/7.66,35.86/) !Parkinson and Washington  
    real, parameter, dimension(2)   :: sens_ice_air = (/0.69e-03,1.12e-03/)       
!	
end module Block_Ice_Snow
