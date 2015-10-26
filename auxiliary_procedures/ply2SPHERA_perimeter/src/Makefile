# Variables to be updated 
VERSION = 8_0_AA
#    bin, debug
EXECUTION = bin
#    -O1 (bin), -g -O0 -fbacktrace -C (gfortran debug),
#    -g -O0 -traceback -C -check bounds -check noarg_temp_created -debug all (ifort debug)
COMPILATION_FLAGS = -O1
#    gfortran, ifort
COMPILER = gfortran
#    -fopenmp, -openmp, (blank for debug)  
OMP_FLAG = 
EXE_DIR = ../$(EXECUTION)/
SOURCES = \
$(wildcard ./*.f90)
# Other variables
MAIN_FILE_ROOT = ./ply2SPHERA_perimeter
CODE = $(EXE_DIR)ply2SPHERA_perimeter_v_$(VERSION)_$(COMPILER)_$(EXECUTION)
OBJECTS = $(patsubst %.f90,%.o,$(SOURCES))
# Targets
#    Correct sequence for compile+link:
#       make touch
#       make
#
#all: compile link remove 
#    "compiling" cannot be launched together with "touching" as ".o" files
#    would not be available yet
#compile: $(SOURCES)
#%.f90: %.o
#	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS) -c 
#link: $(OBJECTS)
#	touch $(CODE).x
#	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS)
#	rm -f $(MAIN_FILE_ROOT).o
all: compile remove
compile: $(SOURCES)
%.f90: %.o
	$(COMPILER) $(MAIN_FILE_ROOT).f90 -o $(CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS)
remove: $(OBJECTS)
	rm -f $^
	rm -f *.mod
touch:
	touch $(OBJECTS)
clean: 
	rm -f $(CODE).x 

