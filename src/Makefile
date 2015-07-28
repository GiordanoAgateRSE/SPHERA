# Variables to be updated 
VERSION = 7_2
DATE = 27Jul15
#    gfortran ifort
COMPILER = gfortran
#    -fopenmp -openmp  
OMP_FLAG = -fopenmp
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
EXE_DIR = ../bin/
CODE = $(EXE_DIR)SPHERA_v_$(VERSION)_$(DATE)
OBJECTS = $(patsubst %.f90,%.o,$(SOURCES))
# Targets
#    Correct sequence for compiling+linking:
#       make touching
#       make 
all: compiling linking removing 
#    "compiling" cannot be launched together with "touching" as ".o" files
#    would not be available yet
compiling: $(SOURCES)
%.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) -O1 -c 
linking: $(OBJECTS)
	touch $(CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(CODE).x $(OMP_FLAG) -O1
	rm -f $(MAIN_FILE_ROOT).o
removing: $(OBJECTS)
	rm -f $^
	rm -f *.mod
touching:
	touch $(OBJECTS)
clean: 
	rm -f $(CODE).x 

