 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
random_huv: $(objects) $(fft_lib) $(sourcedir)/init/random_huv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/random_huv.f90 -o random_huv $(flags)

random_pv: $(objects) $(sourcedir)/init/random_pv.f90
	$(f90) parameters.o constants.o spherical.o randomise.o $(sourcedir)/init/random_pv.f90 -o random_pv $(flags)

random_top: $(objects) $(sourcedir)/init/random_top.f90
	$(f90) parameters.o constants.o spherical.o randomise.o $(sourcedir)/init/random_top.f90 -o random_top $(flags)

bcjet: $(objects) $(sourcedir)/init/bcjet.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/bcjet.f90 -o bcjet $(flags)

axisym_top: $(objects) $(sourcedir)/init/axisym_top.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/axisym_top.f90 -o axisym_top $(flags)

rest: $(objects) $(sourcedir)/init/rest.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest.f90 -o rest $(flags)

ellipole: $(objects) $(sourcedir)/init/ellipole.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/ellipole.f90 -o ellipole $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


