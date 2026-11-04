 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
balinit: $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/init/balinit.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/init/balinit.f90 -o balinit $(flags)

random_top: $(objects) $(fft_lib) $(sourcedir)/init/random_top.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/random_top.f90 -o random_top $(flags)

random_pv: $(objects) $(fft_lib) $(sourcedir)/init/random_pv.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/random_pv.f90 -o random_pv $(flags)

ellipse_top: parameters.o constants.o $(sourcedir)/init/ellipse_top.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/ellipse_top.f90 -o ellipse_top $(flags)

bcjet: parameters.o constants.o $(sourcedir)/init/bcjet.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/bcjet.f90 -o bcjet $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)
