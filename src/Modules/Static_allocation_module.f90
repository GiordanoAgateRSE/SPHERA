!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: Static_allocation_module            
! Description: Module to define global (and statically allocated) variables.                                          
!-------------------------------------------------------------------------------
module Static_allocation_module
! Global constants: start
integer(4),public,parameter :: SPACEDIM = 3
integer(4),public,parameter :: PLANEDIM = 2
! MAXOPENSIDES applies both in 3D and 2D
integer(4),public,parameter :: MAXOPENSIDES = 10
! MAXPARTLINE applies both in 3D and 2D
integer(4),public,parameter :: MAXPARTLINE = 2000
#ifdef SPACE_3D
integer(4),public,parameter :: MAXOPENFACES = 50
integer(4),public,parameter :: MAXFACENODES = 6
integer(4),public,parameter :: INT_KERNELTABLE = 40
#elif defined SPACE_2D
integer(4),public,parameter :: MAXCLOSEBOUNDSIDES = 4
integer(4),public,parameter :: LIMCLOSEBOUNDSIDES = 4
integer(4),public,parameter :: NUMCOLS_BIT = 5
#endif
integer(4),public,parameter :: INIPARTICLEBUFFER = 225000
integer(4),public,parameter :: FI_STEPS = 8
integer(4),public,parameter :: TETA_STEPS = 4 * FI_STEPS
! Global constants for array sizes: end
! Global constants: start
integer(4),public,parameter :: menouno = - 1
double precision,public,parameter :: max_positive_number = 3.4d+38
double precision,public,parameter :: max_negative_number = -3.4d+38
double precision,public,parameter :: PIGRECO = 3.14159265358979d0
double precision,public,parameter :: Nepero_number = 2.718281828d0
! von Karman's constant
double precision,public,parameter :: k_v = 0.41d0
#ifdef SPACE_3D
double precision,public,parameter :: KERNELCONST3D = 1.d0 / PIGRECO
#elif defined SPACE_2D
! 10/(7*PIGRECO)
double precision,public,parameter :: KERNELCONST2D = 0.454728408833987d0 
#endif
double precision,public,parameter :: FI_INTERVAL = PIGRECO / 2.0d0  
double precision,public,parameter :: TETA_INTERVAL = PIGRECO + PIGRECO 
double precision,public,parameter :: zero = 0.d0
double precision,public,parameter :: one = 1.d0 
double precision,public,parameter :: two = 2.d0
double precision,public,parameter :: four = 4.d0 
double precision,public,parameter :: half = 0.5d0 
double precision,public,parameter :: quarter = 0.25d0
! Tolerance for coordinate checking
double precision,public,parameter :: xyz_tolerance= 1.d-4
double precision,public,parameter :: const_m_9999 = -9999.d0 
double precision,public,parameter :: sqrttwo = 1.4142135623731d0
double precision,public,parameter :: arrotondamento = 1.d-5
double precision,public,parameter :: azzeramento = 1.d8
! Gravitational constant
double precision,public,parameter :: GI = 9.80665
! von Karman constant
double precision,public,parameter :: vKconst = 0.41
! Global constants: end
! Global variables: start 
logical :: dt_alfa_Mon
! Flag to pinpoint the second read of the input files
logical :: input_second_read
#ifdef SPACE_3D
! Flag on the presence of CLC input data
logical :: CLC_flag
integer(4),public,parameter :: ncord = 3
#elif defined SPACE_2D
integer(4),public,parameter :: ncord = 2 
#endif
! Maximum file unit booked in the modules
integer(4) :: max_file_unit_booked = 100
integer(4) :: NMedium,NPartZone
integer(4) :: npointst,NPoints,NPointsl,NPointse,NLines
#ifdef SPACE_3D
integer(4) :: GCBFVecDim
#endif
integer(4) :: NumVertici,NumBVertices
! Number of zones
integer(4) :: NumTratti
#ifdef SPACE_3D
integer(4) :: NumFacce
#elif defined SPACE_2D
integer(4) :: NumBSides
#endif
integer(4) :: nag,nagpg,PARTICLEBUFFER
! Variable to count the particles, which are not "sol"
integer(4) :: indarrayFlu
integer(4) :: it_start,on_going_time_step,it_eff,indarraySol
double precision :: simulation_time,elapsed_time,dt,pesodt,dt_average,DTminBER
! Global variables: end
! Global variables for inlet sections: start
! Flag to detect the job-start conditions for the inlet sections
logical :: inlet_job_start_flag = .true.
integer(4) :: mat,irz,izone
#ifdef SPACE_3D
! Number of outlet sections (open sections)
integer(4) :: NumOpenFaces
! ID of the current inlet section or number of inlet sections, depending on 
! the program unit
integer(4) :: SourceFace
#elif defined SPACE_2D
integer(4) :: NumOpenSides,SourceSide
#endif
integer(4) :: itime_jet
double precision :: RowPeriod,yfila,zfila
! Following emission time for inlet particles
double precision :: emission_time
#ifdef SPACE_3D
integer(4),dimension(1:MAXOPENFACES) :: OpenFace
integer(4),dimension(1:MAXOPENSIDES) :: NumPartFace
#elif defined SPACE_2D
integer(4),dimension(1:MAXOPENSIDES) :: NumPartperLine,OpenSide
#endif
double precision,dimension(1:MAXOPENSIDES) :: RowVelocity
double precision,dimension(1:SPACEDIM) :: P
double precision,dimension(1:SPACEDIM) :: Q
double precision,dimension(1:SPACEDIM) :: nn
! This array applies both to 3D and 2D
double precision,dimension(1:MAXOPENSIDES,1:MAXPARTLINE,1:SPACEDIM) :: PartLine
! Global variables for inlet sections: end
! Global variables to compute SA-SPH integrals: start
integer(4),parameter :: BITrows = FI_STEPS * TETA_STEPS
integer(4),parameter :: BITcols = 7
#ifdef SPACE_3D
integer(4),parameter :: ktrows = INT_KERNELTABLE
integer(4),parameter :: ktcols = 4
double precision :: ktdelta
double precision,dimension(0:ktrows,0:ktcols) :: kerneltab
#endif
double precision,dimension(1:BITrows,1:BITcols) :: BoundIntegralTab
! Global variables to compute SA-SPH integrals: end
! Global variables for Paraview output files : start
! Maximum number of blocks
integer(4),public,parameter :: maxnumblock = 9999
! Flag to activate the "vtkconverter"
logical :: vtkconv
! Number of current blocks
integer(4) :: nblocchi,block
! Writing time step
double precision :: freq_time,val_time
! Array to store the number of blocks
integer(4),dimension(maxnumblock) :: blocchi
! Array to store the time of the blocks
double precision,dimension(maxnumblock) :: Time_Block
! Global variables for output files to Paraview: end
integer(4),parameter :: MAXTIT = 10
integer(4),public,parameter :: lencard = 200
! Maximum number of points for the definition of a "GENERIC" area
integer(4),public, parameter :: MAXPOINTSZONE = 20
! Maximum number of data for the definition of a velocity "LAW"
integer(4),public, parameter :: MAXPOINTSVLAW  = 50
! Flag for error file existence in erosion model
logical :: err_flag
! Flag if the run is a restart
logical :: restart
! Flag to kill the execution
logical :: kill_flag
logical :: current_version
#ifdef SOLID_BODIES
! Flags to activate/deactivate pressure limiters on the body surfaces (input)
logical :: body_minimum_pressure_limiter,body_maximum_pressure_limiter
! Flag to activate/deactivate the treatment for submerged boundary thin walls, 
! far from other boundary types
logical :: thin_walls
! Flag to remove/keep fluid particles which have crossed a fluid-body interface
logical :: remove_fluid_in_body
#endif
! Max number of neighbouring particles
integer(4) :: NMAXPARTJ
#ifdef SPACE_3D
! Max number of close boundary faces for the current particle
integer(4) :: MaxNcbf
#elif defined SPACE_2D
! Max number of close boundary sides for the current particle
integer(4) :: MaxNcbs
#endif
#ifdef SOLID_BODIES
! Total number of body particles
integer(4) :: n_body_part
#ifdef SPACE_3D
! Number of CAE-made body particles
integer(4) :: n_body_part_CAE
#endif
! Total number of surface body particles
integer(4) :: n_surf_body_part
! Total number of bodies 
integer(4) :: n_bodies
! Slip conditions for FSI (ref. input file template)
integer(4) :: FSI_slip_conditions
! Number of fluid particles removed from solid bodies after the initial 
! conditions
integer(4) :: fluid_in_body_count
! Number of mass-frozen fluid particles (ALE3)
integer(4) :: mass_frozen_count = 0
#ifdef SPACE_3D
! Number of CAE-made solid bodies
integer(4) :: n_bodies_CAE
#endif
#endif
#ifdef SPACE_3D
! Number of convex edges
integer(4) :: NumBEdges
#endif
! Input variable
double precision :: friction_angle
! 0.001 * Domain%h 
double precision :: eta
! 0.01 * Domain%h * Domain%h
double precision :: eta2
! 2.*Domain%h
double precision :: doubleh
! Domain%h*Domain%h
double precision :: squareh
! doubleh * doubleh
double precision :: square_doubleh
! 1./Domain%h
double precision :: unosuh
! 1./(Domain%h*Domain%h)
double precision :: unosusquareh
#ifdef SOLID_BODIES
! Ratio between fluid particle and body particle size (only for ICs of 
! handmade bodies)
double precision :: dx_dxbodies
! Coefficient for initial removal of fluid particles within solid bodies
double precision :: c_ini_rem_fp_sb
! Numerical times for body dynamics (input)
double precision :: time_max_no_body_gravity_force
double precision :: time_max_no_body_frontier_impingements
#endif
! Indices of cells that must be considered around the current one in the 
! program unit "CalcVarLength"
integer(4),dimension(14,3) :: indicecelle
! Pointer for coordinate location 2D or 3D; = (/0,1,3,0,0,1,2,3/); initialized
! in the main program for compatibility with xlf90 
integer(4),dimension(0:3,2) :: icoordp
character(len=8),parameter :: acode = "SPHERA  "
character(len=8),parameter :: version = "10.0.0"
character(255) :: nomecaso,nomecas2
character(1),dimension(0:3) :: xyzlabel = (/ "T", "X", "Y", "Z" /)  
character(4),dimension(3) :: ncordlabel = (/ "    ", "(2D)", "(3D)" /)
! Killer file name
character(255) :: nomefilekill
! Operating System: "linux" (no other option is active). This string is 
! useful for the real-time assessment of the final elapsed time. Any other 
! string value deactvates the assessment above.
character(10) :: exetype = "linux"
! File name for error file in erosion model
character(255) :: nomefileerr
! "original" or "euristic"
character(len=8) :: dt_opt = "original"
character(len=lencard),dimension(MAXTIT) :: title
end module
