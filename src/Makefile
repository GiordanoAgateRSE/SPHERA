#-------------------------------------------------------------------------------
# SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
# Computational Fluid Dynamics code).
# Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
# formerly CESI-Ricerca di Sistema)
#
# SPHERA authors and email contact are provided in SPHERA documentation.
#
# This file is part of SPHERA.
# SPHERA is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# SPHERA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
# along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# Description: SPHERA Makefile.
#-------------------------------------------------------------------------------

# Variables to be modified: start

#    CODE VERSION
#       Suggested values:
#          " v9_UserInitials_CommitDate"
VERSION = v9_

#    COMPILER
#       Suggested values:
#          " ifort"
#          " gfortran"
COMPILER = gfortran

#    EXECUTION
#       Suggested values:
#          " optimized"
#          " debug"
EXECUTION = debug

#    PARALLELIZATION
#       Suggested values:
#          parallel simulations
#             " OMP"
#          sequential simulations
#             " NO"
PARALLELIZATION = OMP

#    Preprocessor Directives (PD): start

#       SPACE DIMENSIONALITY
#          Suggested values:
#             3D
#                " -DSPACE_3D"
#             2D
#                " -DSPACE_2D"
PD_SPACE = -DSPACE_3D

#       KTGF
#          Suggested values:
#             KTGF scheme
#                " -DKTGF"
#             no KTGF scheme
#                " "
PD_KTGF = -DKTGF

#       SOLID_BODIES
#          Suggested values:
#             Body transport
#                " -DSOLID_BODIES"
#             No body transport
#                " "          
PD_SOLID_BODIES = -DSOLID_BODIES

#    Preprocessor Directives (PD): end

# Variables to be modified: end

# Variable which might be assigned for debug purposes as an alternative to 
# COMPILATION_FLAGS for some program units
COMPILATION_FLAGS_ISSUES =  

# Possible incompatibilities among the functionalities under pre-processor 
# directives
# None

# IDs of the above choices reported in the name of the executable code
ifeq ($(COMPILER), ifort)
   ID_COMPILER = 1
   else ifeq ($(COMPILER), gfortran)
      ID_COMPILER = 2
endif
ifeq ($(EXECUTION), optimized)
   ID_EXECUTION = 1
   else ifeq ($(EXECUTION), debug)
      ID_EXECUTION = 2
endif
ifeq ($(PARALLELIZATION), NO)
   ID_PARALLELIZATION = 0
   else ifeq ($(PARALLELIZATION), OMP)
      ID_PARALLELIZATION = 1
endif
ifeq ($(PD_SPACE), -DSPACE_3D)
   ID_PD_SPACE = 3
   else ifeq ($(PD_SPACE), -DSPACE_2D)
      ID_PD_SPACE = 2
endif
ifeq ($(PD_KTGF), -DKTGF)
   ID_PD_KTGF = 1
   else
      ID_PD_KTGF = 0
endif
ifeq ($(PD_SOLID_BODIES), -DSOLID_BODIES)
   ID_PD_SOLID_BODIES = 1
   else
      ID_PD_SOLID_BODIES = 0
endif

# COMPILATION_FLAGS
COMPILATION_FLAGS = -O2
ifeq ($(EXECUTION), debug)
   ifeq ($(COMPILER), gfortran)
      COMPILATION_FLAGS = -g -O0 -fbacktrace -C -Wall -fbounds-check -ffpe-trap=zero,overflow,underflow,invalid,denormal
      else ifeq ($(COMPILER), ifort)
         COMPILATION_FLAGS = -g -O0 -traceback -C -check bounds -check noarg_temp_created -debug all -fpe0
   endif
endif

# COMPILATION_FLAG_IPO
COMPILATION_FLAG_IPO =  
ifeq ($(EXECUTION), optimized)
   ifeq ($(COMPILER), ifort)
      COMPILATION_FLAG_IPO = -ipo
   endif
endif

# OMP_FLAG
OMP_FLAG =  
ifeq ($(PARALLELIZATION), OMP)
   ifeq ($(COMPILER), ifort)
      OMP_FLAG = -qopenmp
      else ifeq ($(COMPILER), gfortran)
         OMP_FLAG = -fopenmp
   endif
endif

EXE_DIR = ../exe/

SOURCES = \
./Modules/Static_allocation_module.f90 \
./Modules/Hybrid_allocation_module.f90 \
./Modules/Dynamic_allocation_module.f90 \
./Modules/Memory_I_O_interface_module.f90 \
./Modules/Neighbouring_Search_interface_module.f90 \
./Modules/I_O_file_module.f90 \
./Modules/SA_SPH_module.f90 \
./Modules/Time_module.f90 \
$(wildcard ./BC/*.f90) \
$(wildcard ./BE_mass/*.f90) \
$(wildcard ./BE_momentum/*.f90) \
$(wildcard ./Body_Transport/BE_body_momentum/*.f90) \
$(wildcard ./Body_Transport/BE_mass/*.f90) \
$(wildcard ./Body_Transport/BE_momentum/*.f90) \
$(wildcard ./Body_Transport/IC/*.f90) \
$(wildcard ./Body_Transport/Neighbouring_Search/*.f90) \
$(wildcard ./Body_Transport/Post_processing/*.f90) \
$(wildcard ./Body_Transport/Time_Integration/*.f90) \
$(wildcard ./CLC/*.f90) \
$(wildcard ./DB_SPH/*.f90) \
$(wildcard ./EoS/*.f90) \
$(wildcard ./Erosion_Criterion/*.f90) \
$(wildcard ./Geometry/*.f90) \
$(wildcard ./IC/*.f90) \
$(wildcard ./KTGF_tools/*.f90) \
./Main_algorithm/Gest_Trans.f90 \
./Main_algorithm/time_step_loop.f90 \
$(wildcard ./Memory_I_O/*.f90) \
$(wildcard ./Memory_I_O/allocate_de_derived_types_a/*.f90) \
$(wildcard ./Memory_I_O/allocate_de_derived_types_b/*.f90) \
$(wildcard ./Memory_I_O/allocate_de_preset_types/*.f90) \
$(wildcard ./Neighbouring_Search/*.f90) \
$(wildcard ./Post_processing/*.f90) \
$(wildcard ./Post_processing/2D_fields/*.f90) \
$(wildcard ./Post_processing/elapsed_time/*.f90) \
$(wildcard ./Post_processing/monitors/*.f90) \
$(wildcard ./Pre_processing/*.f90) \
$(wildcard ./Pre_processing/reading_inp/*.f90) \
$(wildcard ./SA_SPH/BE_mass/*.f90) \
$(wildcard ./SA_SPH/BE_momentum/*.f90) \
$(wildcard ./SA_SPH/IC/*.f90) \
$(wildcard ./SA_SPH/integrals_IC/*.f90) \
$(wildcard ./SA_SPH/integrals_time_steps/*.f90) \
$(wildcard ./SA_SPH/Neighbouring_Search/*.f90) \
$(wildcard ./SA_SPH/Time_integration/*.f90) \
$(wildcard ./SA_SPH/Turbulence/*.f90) \
$(wildcard ./Strings/*.f90) \
$(wildcard ./Time_Integration/*.f90) 

# Auxiliary variables to be possibly assigned for debug purposes as an 
# alternative to the variable "SOURCES" for some program units
SOURCES_issues =   
SOURCES_b = 
SOURCES_issues_b =
 
MAIN_FILE_ROOT = ./Main_algorithm/sphera

EXE_CODE = \
$(EXE_DIR)SPHERA_$(VERSION)_$(ID_COMPILER)$(ID_EXECUTION)$(ID_PARALLELIZATION)$(ID_PD_SPACE)$(ID_PD_KTGF)$(ID_PD_SOLID_BODIES)

OBJECTS = \
$(patsubst %.f90,%.o,$(SOURCES)$(SOURCES_issues)$(SOURCES_b)$(SOURCES_issues_b))

# Targets: start
#    Correct sequence to compile+link from a Linux terminal: 
#       "make"

all: touch compile compile_issues compile_b compile_issues_b link clean 

compile: $(SOURCES)

compile_b: $(SOURCES_b)

$(SOURCES) $(SOURCES_b): %.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS) $(COMPILATION_FLAG_IPO) -c -cpp $(PD_SPACE) $(PD_KTGF) $(PD_SOLID_BODIES)

compile_issues: $(SOURCES_issues)

compile_issues_b: $(SOURCES_issues_b)

$(SOURCES_issues) $(SOURCES_issues_b): %.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS_ISSUES) -c -cpp $(PD_SPACE) $(PD_KTGF) $(PD_SOLID_BODIES)

link: $(OBJECTS)
	touch $(EXE_CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(EXE_CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS) -cpp $(PD_SPACE) $(PD_KTGF) $(PD_SOLID_BODIES)
	rm -f $(MAIN_FILE_ROOT).o

clean: $(OBJECTS)
	rm -f $^
	rm -f *.mod

touch:
	touch $(OBJECTS)

remove: 
	rm -f $(EXE_CODE).x

# Targets: end
