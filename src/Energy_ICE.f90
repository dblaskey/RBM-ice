      SUBROUTINE ENERGY_ICE   (q_rslt, nd, ncell, nr, ns)
!
Use Block_Energy
Use Block_Hydro
Use Block_Network
Use Block_Ice_Snow
!
Implicit None
integer                   :: ncell,nd,nr,ns
!
real                      :: cnddmmy,cnd_tmp
real                      :: q_ice,q10
real                      :: q_rslt
real                      :: Ltnt_Heat,Sens_Heat
real                      :: cndctvy,delta_ice,LW_back,LW_in,SW_in
real                      :: dvsr,delta_Temp
real                      :: T_B,T_ice,T_p,T_p_cubed,T_riv,T_srfc
real                      :: snow_p,cndctvySI,snow_cndctvy,thickS,T_S
real                      :: delta_snow,lvfs 
!
real,parameter            :: thick0 = 0.8
!
  T_B = 0.0
  T_S = 0.0
  T_srfc = ice_temp(nr,ns,n1)
  T_p = T_srfc + T_Kelvin
  T_p_cubed = T_p*T_p*T_p
!  cndctvy  = ice_cndctvy/(ice_thick(nr,ns,n1)+0.01)
  cndctvy = ice_cndctvy/thick0
  dvsr = 4.0*epsilon*Stfn_Bltz*T_p_cubed + cndctvy 
  LW_back = epsilon*Stfn_Bltz*T_p*T_p_cubed 
  LW_in = epsilon*QNA(ncell)

!
  call Vapor(q_ice,q10,T_p,ncell,ns)
  Ltnt_Heat = rho_air*lvs*Ch_Cg*wind(ncell)*(q10-q_ice)
!
  Sens_Heat = rho_air*cp_air*Ch_Cg*wind(ncell)*(dbt(ncell)-T_srfc)
!

if (snow_p .gt. 0.01) then ! Need to allocate snow_p
  SW_in = (1.0-snow_albedo)*QNS(ncell)
  cndctvySI = ice_cndctvy*snow_cndctvy/(thick0*snow_cndctvy)               & ! Allocated cndctvySI, snow_cndctvy, thickS
              + thickS*ice_cndctvy
  dvsr = 4.0*epsilon*Stfn_Bltz*T_p_cubed + cndctvySI 
  delta_Temp = (Sens_Heat + Ltnt_Heat + LW_in + SW_in - LW_back           &
            +cndctvySI*(T_B-T_srfc))/dvsr
!  
  snow_temp(nr,ns,n2) = snow_temp(nr,ns,n1) + delta_Temp !ALLOCATE SNOW TEMP
  snow_temp(nr,ns,n2) = AMAX1(-.01,snow_temp(nr,ns,n2))
  ! Check here to see if snow is thawing 
  !
  if (snow_temp(nr,ns,n2) .ge. T_S) then !ALLOCATE T_S
    ! NEED TO DOUBLE CHECK THIS IS RIGHT
    delta_snow = dt_comp*(Sens_Heat + Ltnt_Heat + LW_in + SW_in - LW_back)/lvfs ! ALLOCATE DELTA SNOW, LVFS
    snow_thick(nr,ns,n2) = snow_thick(nr,ns,n1) - delta_snow                
    snow_thick(nr,ns,n2) = AMAX1(snow_thick(nr,ns,n2),0.0)
  else ! SHOULD AT LEAST UPDATE TO ACCOUNT FOR INSULATION OF ICE (changed to conductivity of snow and ice)
    T_ice = ice_temp(nr,ns,n1)
    T_riv = temp(nr,ns,n1)
    delta_ice = dt_comp*cndctvySI*(T_riv-T_ice)/lvf
    ice_thick(nr,ns,n2) = ice_thick(nr,ns,n1) + delta_ice
    ice_thick(nr,ns,n2) = AMIN1(ice_thick(nr,ns,n2),depth(ncell))
    ICE(ncell) = 200.
  end if
else
  SW_in = (1.0-0.4*I0)*(1.0-ice_albedo)*QNS(ncell) ! Check if this is NEEDED in Vapor, may need to move all of Vapor
  
  delta_Temp = (Sens_Heat + Ltnt_Heat + LW_in + SW_in - LW_back           &
            +cndctvy*(T_B-T_srfc))/dvsr
  cnd_tmp = cndctvy*(T_B-T_srfc)


  ice_temp(nr,ns,n2) = ice_temp(nr,ns,n1) + delta_Temp
  ice_temp(nr,ns,n2) = AMAX1(-.01,ice_temp(nr,ns,n2))
!
!
! Check here to see if ice is thawing UPDATE TO MOVE THIS INSIDE IF ELSE SNOW STATEMENT
!
  if (ice_temp(nr,ns,n2) .ge. T_B) then
    
    delta_ice = dt_comp*(Sens_Heat + Ltnt_Heat + LW_in + SW_in - LW_back)/lvf
   ice_thick(nr,ns,n2) = ice_thick(nr,ns,n1) - delta_ice                
    ice_thick(nr,ns,n2) = AMAX1(ice_thick(nr,ns,n2),0.0)
    ice_temp(nr,ns,n2) = AMIN1(ice_temp(nr,ns,n2),T_B)
    ice_temp(nr,1:ns,n2) = T_B
    ICE(ncell) = 102.
!
! Freezing
!
  else
    T_ice = ice_temp(nr,ns,n1)
    T_riv = temp(nr,ns,n1)
    delta_ice = dt_comp*cndctvy*(T_riv-T_ice)/lvf
    ice_thick(nr,ns,n2) = ice_thick(nr,ns,n1) + delta_ice
    ice_thick(nr,ns,n2) = AMIN1(ice_thick(nr,ns,n2),depth(ncell))
    ICE(ncell) = 200.
  end if
end if
!
! If ice thickness is less than the minimum, river is no longer frozen
!
    if (ice_thick(nr,ns,n2) .lt. 0.01) then 
      ICE(ncell) = 101.
      ice_thick(nr,ns,n2) = 0.00
    end if
!
else
  
!
! If ice thickness is less than the minimum, river is no longer frozen
!
    if (ice_thick(nr,ns,n2) .lt. 0.01) then 
      ICE(ncell) = 101.
      ice_thick(nr,ns,n2) = 0.00
    end if
!
!
  if (ice_thick(nr,ns,n2) .gt. depth(ncell)) ice_thick(nr,ns,n2) = depth(ncell)
!
  q_rslt =  SW_in + LW_in -LW_back - Ltnt_Heat + Sens_Heat

  
      RETURN
      END
