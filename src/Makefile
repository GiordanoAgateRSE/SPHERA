# Variables to be updated 
VERSION = 8_0_AA
#    bin, debug
EXECUTION = bin
#    -O1 (bin), -g -O0 -fbacktrace -C (gfortran debug),
#    -g -O0 -traceback -C -check bounds -check noarg_temp_created -debug all (ifort debug)
COMPILATION_FLAGS = -O1
#    gfortran, ifort
COMPILER = ifort
#    -fopenmp, -openmp, (blank for debug)  
OMP_FLAG = -openmp
EXE_DIR = ../$(EXECUTION)/
SOURCES = \
./Modules/Static_allocation_module.f90 \
./Modules/Hybrid_allocation_module.f90 \
./Modules/Dynamic_allocation_module.f90 \
./Modules/I_O_diagnostic_module.f90 \
./Modules/I_O_ENG_module.f90 \
./Modules/I_O_file_module.f90 \
./Modules/I_O_ITA_module.f90 \
./Modules/I_O_language_module.f90 \
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
./Main_algorithm/Loop_Irre_2D.f90 \
./Main_algorithm/Loop_Irre_3D.f90 \
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
#    Correct sequence for compile+link:
#       make touch
#       make 
all: compile link remove 
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

