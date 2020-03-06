# 6 compilation procedures have to be tested (without warnings):
#    gfortran_debug_omp, gfortran_debug, gfortran_bin, ifort_debug_omp, 
#    ifort_debug, ifort_bin 
# Variables to be updated
VERSION = 9_0_0
# Suggested values for EXECUTION
#    bin
#    debug
#    debug_omp
EXECUTION = debug_omp
# Suggested values for COMPILATION_FLAGS
#   -O1 (ifort_bin and gfortran_bin)
#   -g -O0 -fbacktrace -C -Wall -fbounds-check (gfortran debug and gfortran_debug_omp),
#   -g -O0 -traceback -C -check bounds -check noarg_temp_created -debug all -fpe0 (ifort_debug and ifort_debug_omp)
COMPILATION_FLAGS = -g -O0 -fbacktrace -C -Wall -fbounds-check
# Suggested values for COMPILER
#    gfortran
#    ifort
COMPILER = gfortran
# Suggested values for OMP_FLAG
#    -fopenmp (gfortran_bin and gfortran_debug_omp)
#    -qopenmp (ifort_bin and ifort_debug_omp) 
OMP_FLAG = -fopenmp
EXE_DIR = ../$(EXECUTION)/
SOURCES = \
./Modules/Static_allocation_module.f90 \
./Modules/Hybrid_allocation_module.f90 \
./Modules/Dynamic_allocation_module.f90 \
./Modules/I_O_diagnostic_module.f90 \
./Modules/I_O_file_module.f90 \
./Modules/SA_SPH_module.f90 \
./Modules/Time_module.f90 \
$(wildcard ./BC/*.f90) \
$(wildcard ./BE_mass/*.f90) \
$(wildcard ./BE_momentum/*.f90) \
$(wildcard ./Body_Transport/*.f90) \
$(wildcard ./Constitutive_Equation/*.f90) \
$(wildcard ./DB_SPH/*.f90) \
$(wildcard ./Erosion_Criterion/*.f90) \
$(wildcard ./Geometry/*.f90) \
$(wildcard ./IC/*.f90) \
$(wildcard ./Interface_dispersion/*.f90) \
./Main_algorithm/Gest_Dealloc.f90 \
./Main_algorithm/Gest_Trans.f90 \
./Main_algorithm/time_step_loop.f90 \
$(wildcard ./Neighbouring_Search/*.f90) \
$(wildcard ./Post_processing/*.f90) \
$(wildcard ./Pre_processing/*.f90) \
$(wildcard ./SA_SPH/*.f90) \
$(wildcard ./Strings/*.f90) \
$(wildcard ./Time_Integration/*.f90)
# Other variables
MAIN_FILE_ROOT = ./Main_algorithm/sphera
CODE = $(EXE_DIR)SPHERA_v_$(VERSION)_$(COMPILER)_$(EXECUTION)
OBJECTS = $(patsubst %.f90,%.o,$(SOURCES))
# Targets
#    Correct sequence to compile+link:
#       make touch
#       make
all: touch compile link remove 
#    "compiling" cannot be launched together with "touching" as ".o" files
#    would not be available yet
compile: $(SOURCES)
%.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS) -c 
link: $(OBJECTS)
	touch $(CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS)
	rm -f $(MAIN_FILE_ROOT).o
remove: $(OBJECTS)
	rm -f $^
	rm -f *.mod
touch:
	touch $(OBJECTS)
clean: 
	rm -f $(CODE).x 

