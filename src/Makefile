#-------------------------------------------------------------------------------
# SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
# Computational Fluid Dynamics code).
# Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
#          " v9_UserInitials_CommitData"
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

#    Pre-Processor Directives (PPD): start

#       DIMENSIONALITY
#          Suggested values:
#             3D
#                " -DDIMENSIONALITY=3D"
#             2D
#                " -DDIMENSIONALITY=2D"
PPD_DIMENSIONALITY = -DDIMENSIONALITY=2D

#       VERBOSITY
#          Suggested values:
#             Synthetic output
#                " "
#             Verbose output
#                " -DVERBOSITY"
PPD_VERBOSITY =  

#       KTGF
#          Suggested values:
#             KTGF scheme with possible 2-interface 3D erosion criterion
#                " -DKTGF=FULL"
#             no KTGF scheme
#                " -DKTGF=NO"
#             KTGF scheme with possible 1-interface 2D erosion criterion
#                " -DKTGF=EC2D"
PPD_KTGF = -DKTGF=NO

#       SOLID_BODIES
#          Suggested values:
#             Body transport
#                " -DSOLID_BODIES"
#             No body transport
#                " "          
PPD_SOLID_BODIES = -DSOLID_BODIES

#       BOUNDARY_TREATMENT
#          Suggested values:
#             SASPH
#                " -DBOUNDARY_TREATMENT=SASPH"
#             DBSPH
#                " -DBOUNDARY_TREATMENT=DBSPH"
PPD_BOUNDARY_TREATMENT = -DBOUNDARY_TREATMENT=SASPH

#       TIME_INTEGRATION
#          Suggested values:
#             Leapfrog scheme
#                " -DTIME_INTEGRATION=LEAPFROG"
#             Euler scheme
#                " -DTIME_INTEGRATION=EULER"
#             Heun scheme
#                " -DTIME_INTEGRATION=HEUN"
PPD_TIME_INTEGRATION = -DTIME_INTEGRATION=LEAPFROG

#    Pre-Processor Directives (PPD): end

# Variables to be modified: end

# Variable which might be assigned for debug purposes as an alternative to 
# COMPILATION_FLAGS for some program units
COMPILATION_FLAGS_issues =  

# Code scheme incompatibilities (PPD)
ifneq ($(PPD_KTGF), -DKTGF=NO)
   PPD_BOUNDARY_TREATMENT = -DBOUNDARY_TREATMENT=SASPH
endif
ifeq ($(PPD_SOLID_BODIES), -DSOLID_BODIES)
   PPD_BOUNDARY_TREATMENT = -DBOUNDARY_TREATMENT=SASPH
   PPD_TIME_INTEGRATION = -DTIME_INTEGRATION=LEAPFROG
endif
ifeq ($(PPD_BOUNDARY_TREATMENT), -DBOUNDARY_TREATMENT=DBSPH)
   PPD_TIME_INTEGRATION = -DTIME_INTEGRATION=EULER
endif

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
ifeq ($(PPD_DIMENSIONALITY), -DDIMENSIONALITY=3D)
   ID_PPD_DIMENSIONALITY = 3
   else ifeq ($(PPD_DIMENSIONALITY), -DDIMENSIONALITY=2D)
      ID_PPD_DIMENSIONALITY = 2
endif
ifeq ($(PPD_VERBOSITY), -DVERBOSITY)
   ID_PPD_VERBOSITY = 1
   else
      ID_PPD_VERBOSITY = 0
endif
ifeq ($(PPD_KTGF), -DKTGF=FULL)
   ID_PPD_KTGF = 1
   else ifeq ($(PPD_KTGF), -DKTGF=NO)
      ID_PPD_KTGF = 0
      else ifeq ($(PPD_KTGF), -DKTGF=EC2D)
         ID_PPD_DIMENSIONALITY = 2
endif
ifeq ($(PPD_SOLID_BODIES), -DSOLID_BODIES)
   ID_PPD_SOLID_BODIES = 1
   else
      ID_PPD_SOLID_BODIES = 0
endif
ifeq ($(PPD_BOUNDARY_TREATMENT), -DBOUNDARY_TREATMENT=SASPH)
   ID_PPD_BOUNDARY_TREATMENT = 1
   else ifeq ($(PPD_BOUNDARY_TREATMENT), -DBOUNDARY_TREATMENT=DBSPH)
      ID_PPD_DIMENSIONALITY = 2
endif
ifeq ($(PPD_TIME_INTEGRATION), -DTIME_INTEGRATION=EULER)
   ID_PPD_TIME_INTEGRATION = 1
   else ifeq ($(PPD_TIME_INTEGRATION), -DTIME_INTEGRATION=HEUN)
      ID_PPD_TIME_INTEGRATION = 2
      else ifeq ($(PPD_TIME_INTEGRATION), -DTIME_INTEGRATION=LEAPFROG)
         ID_PPD_TIME_INTEGRATION = 3
endif

# COMPILATION_FLAGS
COMPILATION_FLAGS = -O2
ifeq ($(EXECUTION), debug)
   ifeq ($(COMPILER), gfortran)
      COMPILATION_FLAGS = -g -O0 -fbacktrace -C -Wall -fbounds-check
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
./Modules/I_O_diagnostic_module.f90 \
./Modules/I_O_file_module.f90 \
./Modules/SA_SPH_module.f90 \
./Modules/Time_module.f90 \
$(wildcard ./BC/*.f90) \
$(wildcard ./BE_mass/*.f90) \
$(wildcard ./BE_momentum/*.f90) \
$(wildcard ./Body_Transport/*.f90) \
$(wildcard ./DB_SPH/*.f90) \
$(wildcard ./Erosion_Criterion/*.f90) \
$(wildcard ./Geometry/*.f90) \
$(wildcard ./IC/*.f90) \
$(wildcard ./KTGF_tools/*.f90) \
./Main_algorithm/Gest_Dealloc.f90 \
./Main_algorithm/Gest_Trans.f90 \
./Main_algorithm/time_step_loop.f90 \
$(wildcard ./Neighbouring_Search/*.f90) \
$(wildcard ./Post_processing/*.f90) \
$(wildcard ./Pre_processing/*.f90) \
$(wildcard ./SA_SPH/BE_mass/*.f90) \
$(wildcard ./SA_SPH/BE_momentum/*.f90) \
$(wildcard ./SA_SPH/IC/*.f90) \
$(wildcard ./SA_SPH/integrals_IC/*.f90) \
$(wildcard ./SA_SPH/integrals_time_steps/*.f90) \
$(wildcard ./SA_SPH/Neighbouring_Search/*.f90) \
$(wildcard ./SA_SPH/Time_integration/*.f90) \
$(wildcard ./Strings/*.f90) \
$(wildcard ./Time_Integration/*.f90) 

# Auxiliary variables to be possibly assigned for debug purposes as an 
# alternative to the variable "SOURCES" for some program units
SOURCES_issues =   
SOURCES_b = 
SOURCES_issues_b =
 
MAIN_FILE_ROOT = ./Main_algorithm/sphera

EXE_CODE = \
$(EXE_DIR)SPHERA_$(VERSION)_$(ID_COMPILER)$(ID_EXECUTION)$(ID_PARALLELIZATION)$(ID_PPD_DIMENSIONALITY)$(ID_PPD_VERBOSITY)$(ID_PPD_KTGF)$(ID_PPD_SOLID_BODIES)$(ID_PPD_BOUNDARY_TREATMENT)$(ID_PPD_TIME_INTEGRATION)

OBJECTS = \
$(patsubst %.f90,%.o,$(SOURCES)$(SOURCES_issues)$(SOURCES_b)$(SOURCES_issues_b))

# Targets: start
#    Correct sequence to compile+link from a Linux terminal: 
#       "make"

all: touch compile compile_issues compile_b compile_issues_b link remove 

compile: $(SOURCES)

compile_b: $(SOURCES_b)

$(SOURCES) $(SOURCES_b): %.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS_O2) $(COMPILATION_FLAG_IPO) -c $(PPD_DIMENSIONALITY) $(PPD_VERBOSITY) $(PPD_KTGF) $(PPD_SOLID_BODIES) $(PPD_BOUNDARY_TREATMENT) $(PPD_TIME_INTEGRATION)

compile_issues: $(SOURCES_issues)

compile_issues_b: $(SOURCES_issues_b)

$(SOURCES_issues) $(SOURCES_issues_b): %.f90: %.o
	$(COMPILER) $@ -o $< $(OMP_FLAG) $(COMPILATION_FLAGS_O1) -c $(PPD_DIMENSIONALITY) $(PPD_VERBOSITY) $(PPD_KTGF) $(PPD_SOLID_BODIES) $(PPD_BOUNDARY_TREATMENT) $(PPD_TIME_INTEGRATION)

link: $(OBJECTS)
	touch $(EXE_CODE).x
	$(COMPILER) $(MAIN_FILE_ROOT).f90 $^ -o $(EXE_CODE).x $(OMP_FLAG) $(COMPILATION_FLAGS_O2) $(PPD_DIMENSIONALITY) $(PPD_VERBOSITY) $(PPD_KTGF) $(PPD_SOLID_BODIES) $(PPD_BOUNDARY_TREATMENT) $(PPD_TIME_INTEGRATION)
	rm -f $(MAIN_FILE_ROOT).o

remove: $(OBJECTS)
	rm -f $^
	rm -f *.mod

touch:
	touch $(OBJECTS)

clean: 
	rm -f $(EXE_CODE).x

# Targets: end
