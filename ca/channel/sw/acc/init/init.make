 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
geobal: $(objects) $(fft_lib) $(sourcedir)/init/geobal.f90
	$(f90) $(fft_lib) parameters.o constants.o generic.o spectral.o $(sourcedir)/init/geobal.f90 -o geobal $(flags)

igw: $(objects) $(sourcedir)/init/igw.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/igw.f90 -o igw $(flags)

bump: $(objects) $(sourcedir)/init/bump.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/bump.f90 -o bump $(flags)

# Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


