 # Set existence of directory variable used in main makefile:
init_exists = true
 # Calculate f90 codes existing in init directory for making
 # with 'make all': 
present_init_files = $(notdir $(basename $(wildcard $(sourcedir)/init/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:

inigamma: $(objects) $(fft_lib) $(sourcedir)/init/inigamma.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/init/inigamma.f90 -o inigamma $(flags)

ranpv: $(objects) $(sourcedir)/init/ranpv.f90
	$(f90) parameters.o constants.o spherical.o randomise.o $(sourcedir)/init/ranpv.f90 -o ranpv $(flags)

ranhh: $(objects) $(sourcedir)/init/ranhh.f90
	$(f90) parameters.o constants.o spherical.o randomise.o $(sourcedir)/init/ranhh.f90 -o ranhh $(flags)

topo: $(objects) $(sourcedir)/init/topo.f90
	$(f90) parameters.o constants.o spherical.o randomise.o $(sourcedir)/init/topo.f90 -o topo $(flags)

rest: $(objects) $(sourcedir)/init/rest.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/rest.f90 -o rest $(flags)

dambreak: $(objects) $(sourcedir)/init/dambreak.f90
	$(f90) parameters.o constants.o $(sourcedir)/init/dambreak.f90 -o dambreak $(flags)

rossby: $(sourcedir)/init/rossby.f90
	$(f90) $(sourcedir)/init/rossby.f90 -o rossby $(flags)

 # Phony definitions:
.PHONY: init_all
 # Rule for 'make all' in the main make file:
init_all: $(present_init_files)


