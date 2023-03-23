module Block_Energy
!
!   Averaging period (days)
  integer, parameter                  :: ishft = -1,nn_avg = 7
!
!   Energy budget variables
!
!   Incoming short wave radiation, Watts/meter**2
!
  real, dimension(:), allocatable     :: qns
!
!   Incoming atmospheric radiation, Watts/meter**2
!
  real, dimension(:), allocatable     :: qna
!
!   Air temperature at surface, deg. C
!
  real, dimension(:), allocatable     :: dbt
!  
!   Wind speed, m/sec
!
  real, dimension(:), allocatable     :: wind
!
!   Vapor pressure of air at surface, kPa
!
  real, dimension(:), allocatable     :: ea
!
!   Air pressure at surface, kPa
!  
  real, dimension(:), allocatable     :: press 
!
!   Snow depth, m
!  
  real, dimension(:), allocatable     :: snow_p 
!
! Simulated water temperatures 
!   
  real, dimension(:),   allocatable   :: T_head,T_trib 
  real, dimension(:),   allocatable   :: tmp_arry 
  real, dimension(:,:), allocatable   :: Q_nps
  real, dimension(:,:), allocatable   :: T_smth,temp_trib
  real, dimension(:,:), allocatable   :: temp_nps,thermal
    !
  real, dimension(:,:,:), allocatable :: temp

!
  real, dimension (:), allocatable   :: mu,alphamu,beta,gmma,smooth_param
!
  real, dimension (:,:), allocatable :: f1_tilde
!   Some important constants
!
  real             :: kcal_Wsec = 4184.0
  real             :: dlta,dlta1,dlta2,T_0,t1,t2,x1,x2
  real,parameter   :: alpha_ice = 0.5        ! Ice reflectivity - Parkinson-Washington
  real,parameter   :: epsilon = 0.97         ! Emissivity
  real,parameter   :: Stfn_Bltz = 5.67e-08   ! Stefan-Boltzmann - W/m**2/degK
  real,parameter   :: cp_air = 1004.0        ! J/Kg/deg K 
  real,parameter   :: kappa_ice = 8.0        ! J/sec/m**2/DegK - Wanders et al
  real,parameter   :: lvp = 2.501e06         ! Latent heat of vaorization - J(W sec)/kg
  real,parameter   :: lvf = 0.334e06         ! Latent heat of fusion - J(W sec)/kg
  real,parameter   :: lvs = 2.834e06         ! Latent heat of sublimation                    
  real,parameter   :: comp_eps = 0.378       ! Complement of eps_air_h2o
  real,parameter   :: eps_air_h2o = 0.622    ! Ratio 
  real,parameter   :: rho = 1000.0           ! Kg/meter**3
  real,parameter   :: wind_fctr = 1.0        ! dimensionless
  real,parameter   :: rho_Cp = 1000.0*4184.0 ! (Kg/m**3)*(J/(Kg*degK)
  real,parameter   :: rho_air = 1.292        ! Kg/m**3
  real,dimension(2),parameter :: evrate = (/1.5e-8,1.5e-8/) ! 1/kPa
  real,dimension(13),parameter :: alpha_h2o = (/0.457,0.267,0.144,0.091,0.074,0.070,0.072  &
                                             ,0.084,0.124,0.227,0.417,0.542,0.542/) !monthly h20 albedo
!
end module Block_Energy  
