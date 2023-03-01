#
# Makefile for the grid-based semi-Lagrangian water 
# temperature model, RBM10_VIC, modified to predict freezing and thawing
#
# Start of the makefile
#
# Defining variables
#
objects = rbm10_mizu.o Begin.o Systmm.o Particle_Track.o \
	Julian.o Leap_Year.o Read_Forcing.o \
	Ice_Free.o Energy_ICE.o Energy_NO_ICE.o \
	Block_Energy.o Block_Hydro.o Block_Ice_Snow.o Block_Network.o\
	Vapor.o Water_Balance.o Write.o
# GNU F90 compiler	  
f90comp = gfortran -Wall -g3 -fbounds-check
# Makefile
rbm10_VIC: $(objects)
	$(f90comp) -o rbm10_mizu $(objects)
Block_Energy.o: Block_Energy.f90
	$(f90comp) -c Block_Energy.f90
block_energy.mod: Block_Energy.f90
	$(f90comp) -c Block_Energy.f90
Block_Hydro.o: Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90
block_hydro.mod: Block_Hydro.o Block_Hydro.f90
	$(f90comp) -c Block_Hydro.f90
Block_Ice_Snow.o: Block_Ice_Snow.f90
	$(f90comp) -c Block_Ice_Snow.f90
block_ice_snow.mod: Block_Ice_Snow.o Block_Ice_Snow.f90
	$(f90comp) -c Block_Ice_Snow.f90
Block_Network.o: Block_Network.f90
	$(f90comp) -c Block_Network.f90
block_network.mod: Block_Network.o Block_Network.f90
	$(f90comp) -c Block_Network.f90
Begin.o: block_energy.mod block_network.mod block_hydro.mod \
	block_ice_snow.mod Begin.f90
	$(f90comp) -c Begin.f90
Read_Forcing.o: block_hydro.mod block_network.mod Read_Forcing.f90
	$(f90comp) -c Read_Forcing.f90 
Systmm.o: block_network.mod block_hydro.mod block_energy.mod\
	block_ice_snow.mod Systmm.f90
	$(f90comp) -c Systmm.f90
Ice_Free.o: block_hydro.mod  block_network.mod block_energy.mod \
	block_ice_snow.mod Ice_Free.f90
	$(f90comp) -c Ice_Free.f90
Energy_ICE.o: block_hydro.mod  block_network.mod block_energy.mod \
	block_ice_snow.mod Energy_ICE.f90
	$(f90comp) -c Energy_ICE.f90
Energy_NO_ICE.o: block_hydro.mod  block_network.mod block_energy.mod \
	block_ice_snow.mod Energy_NO_ICE.f90
	$(f90comp) -c Energy_NO_ICE.f90
Vapor.o: block_energy.mod block_ice_snow.mod Vapor.f90
	$(f90comp) -c Vapor.f90
Particle_Track.o: block_hydro.mod block_network.mod Particle_Track.f90
	$(f90comp) -c Particle_Track.f90
Write.o: Write.f90 block_ice_snow.mod block_network.mod
	$(f90comp) -c Write.f90
Julian.o: Julian.f90
	$(f90comp) -c Julian.f90
Leap_Year.o: Leap_Year.f90
	$(f90comp) -c Leap_Year.f90
#tntrp.o: tntrp.f90#
#	$(f90comp) -c tntrp.f90
Water_Balance.o: block_hydro.mod block_network.mod Water_Balance.f90
	$(f90comp) -c Water_Balance.f90
rbm10_VIC.o: rbm10_VIC.f90
	$(f90comp) -c rbm10_mizu.f90

# Cleaning everything
clean:
	rm block_energy.mod block_hydro.mod block_network.mod\
           rbm10_mizu
	rm $(objects)
