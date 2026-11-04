 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
fopvbal: $(objects) $(fft_lib) $(sourcedir)/post/fopvbal.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/fopvbal.f90 -o fopvbal $(flags)

froude: $(objects) $(fft_lib) $(sourcedir)/post/froude.f90
	$(f90) $(fft_lib) parameters.o constants.o spectral.o $(sourcedir)/post/froude.f90 -o froude $(flags)

# Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)
