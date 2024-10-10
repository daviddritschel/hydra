 # Set existence of directory variable used in main makefile:
post_exists = true
 # Calculate f90 codes existing in post directory for making
 # with 'make all': 
present_post_files = $(notdir $(basename $(wildcard $(sourcedir)/post/*.f90)))

#---------------------------------------------------------------------------------
 #Rules:
genfg: $(objects) $(fft_lib) $(sourcedir)/post/genfg.f90
	$(f90) $(fft_lib) $(objects) $(sourcedir)/post/genfg.f90 -o genfg $(flags)

 # Phony definitions:
.PHONY: post_all
 # Rule for 'make all' in the main make file:
post_all: $(present_post_files)
