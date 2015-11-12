!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: Static_allocation_module            
! Description: Module to define global (and statically allocated) variables.                                          
!----------------------------------------------------------------------------------------------------------------------------------

module Static_allocation_module
! Global constants for array sizes: start
integer(4),public,parameter :: SPACEDIM = 3
integer(4),public,parameter :: PLANEDIM = 2
integer(4),public,parameter :: MAXCLOSEBOUNDSIDES = 2 
integer(4),public,parameter :: MAXOPENSIDES = 10
integer(4),public,parameter :: MAXOPENFACES = 50
integer(4),public,parameter :: MAXPARTLINE = 2000
integer(4),public,parameter :: LIMCLOSEBOUNDSIDES = 2
integer(4),public,parameter :: MAXFACENODES = 4
integer(4),public,parameter :: NUMCOLS_BIT = 5
integer(4),public,parameter :: NUMROWS_BIT = 22
integer(4),public,parameter :: INT_KERNELTABLE = 40
integer(4),public,parameter :: INIPARTICLEBUFFER = 225000
integer(4),public,parameter :: FI_STEPS = 8
integer(4),public,parameter :: TETA_STEPS = 4 * FI_STEPS
! Global constants for array sizes: end
! Global constants: start
integer(4),public,parameter :: menouno = - 1
double precision,public,parameter :: max_positive_number = 3.4d+38
double precision,public,parameter :: max_negative_number = -3.4d+38
double precision,public,parameter :: PIGRECO = 3.14159265358979 ! PI GRECO
! =10/(7*PIGRECO) 
double precision,public,parameter :: KERNELCONST2D = 0.454728408833987d0 
double precision,public,parameter :: KERNELCONST3D = 1.0d0 / PIGRECO 
double precision,public,parameter :: FI_INTERVAL = PIGRECO/2.0d0  
double precision,public,parameter :: TETA_INTERVAL = PIGRECO + PIGRECO 
double precision,public,parameter :: zero = 0.0d0  
double precision,public,parameter :: one = 1.0d0  
double precision,public,parameter :: two = 2.0d0  
double precision,public,parameter :: four = 4.0d0 
double precision,public,parameter :: half = 0.5d0 
double precision,public,parameter :: quarter = 0.25d0 
double precision,public,parameter :: xyz_tolerance= 1.0d-4 ! tolerance for 
                                                           !coordinate checking
double precision,public,parameter :: const_m_9999 = -9999.0d0 
double precision,public,parameter :: sqrttwo = 1.4142135623731d0 ! dsqrt(2.d0)
double precision,public,parameter :: arrotondamento = 1.0d-5 ! tolerance
double precision,public,parameter :: azzeramento = 1.0d8 ! tolerance
double precision,public,parameter :: GI = 9.80665 ! gravitational constant
double precision,public,parameter :: vKconst = 0.41 ! von Karman constant
! Global constants: end
! Global variables: start 
logical :: dt_alfa_Mon
integer(4) :: Ncord,NMedium,NPartZone
integer(4) :: NPointst,NPoints,NPointsl,NPointse,NLines,NSections,GCBFVecDim
integer(4) :: NumVertici,NumFacce,NumTratti,NumBVertices,NumBSides
integer(4) :: Nag,nagpg,PARTICLEBUFFER,nag_reservoir_CartTopog
integer(4) :: it_start,it_corrente,it_eff,indarrayFlu,indarraySol
double precision :: tempo,dt,pesodt,dt_average,DTminBER
! Global variables: end
! Global variables for inlet sections: start
integer(4) :: NumOpenSides,SourceSide,mat,irz,izone
integer(4) :: SourceFace,NumOpenFaces,pinttimeratio,itime_jet
double precision :: RowPeriod,yfila,ParticleVolume,zfila
integer(4),dimension(1:MAXOPENSIDES) :: NumPartperLine, NumPartFace
integer(4),dimension(1:MAXOPENSIDES) :: Openside
integer(4),dimension(1:MAXOPENFACES) :: OpenFace
double precision,dimension(1:MAXOPENSIDES) :: RowVelocity
double precision,dimension(1:SPACEDIM) :: P
double precision,dimension(1:SPACEDIM) :: Q
double precision,dimension(1:SPACEDIM) :: nn
double precision,dimension(1:MAXOPENSIDES,1:MAXPARTLINE,1:SPACEDIM) :: PartLine
! Global variables for inlet sections: end
! Global variables to compute SA-SPH integrals: start
integer(4),parameter :: BITrows = FI_STEPS * TETA_STEPS
integer(4),parameter :: BITcols = 7
integer(4),parameter :: ktrows = INT_KERNELTABLE
integer(4),parameter :: ktcols = 4
double precision :: ktdelta
double precision,dimension(0:ktrows,0:ktcols) :: kerneltab
double precision,dimension(1:BITrows,1:BITcols) :: BoundIntegralTab
! Global variables to compute SA-SPH integrals: end
! Global variables for output files to Paraview: start
integer(4),public,parameter :: maxnumblock = 9999 ! Maximum number of blocks
logical :: vtkconv ! Flag to activate the "vtkconverter"
integer(4) :: nblocchi,block ! Number of current blocks 
double precision :: freq_time, val_time ! Writing time step 
integer(4),dimension(maxnumblock) :: blocchi ! Array to store the number of 
                                             ! blocks
double precision,dimension(maxnumblock) :: Time_Block ! Array to store the
                                                      ! time of the blocks
! Global variables for output files to Paraview: end
! Global variables for the diffusion model: start 
double precision,public,parameter :: VFmn = 0.1d0 ! Minimum volume fraction 
double precision,public,parameter :: VFmx = 0.9d0 ! Maximum volume fraction
! Global variables for the diffusion model: end
integer(4),parameter :: MAXTIT = 10
integer(4),public,parameter :: lencard= 200
integer(4),public, parameter :: MAXPOINTSZONE = 20 ! Maximum number of points 
                                                   ! for the definition of a
                                                   ! "GENERIC" area
integer(4),public, parameter :: MAXPOINTSVLAW  = 50 ! Maximum number of data 
                                                    ! for the definition of 
                                                    ! a velocity "LAW"
logical :: err_flag ! Flag for error file existence in erosion model
logical :: Restart ! Flag if the run is a restart
logical :: kill_flag ! Flag to kill the execution
logical :: current_version   
logical :: diffusione ! Flag to activate the diffusion model
logical :: esplosione ! Flag to activate the explosion model
logical :: erosione ! Flag to activate the erosion criterion 
integer(4) :: NMAXPARTJ ! Max number of particles surrounding the current one
integer(4) :: MaxNcbs ! Max number of close boundary sides for the current 
                      ! particle
integer(4) :: MaxNcbf ! Max number of close boundary faces for the current 
                      ! particle
integer(4) :: n_body_part ! Total number of body particles
integer(4) :: n_surf_body_part ! Total number of surface body particles 
integer(4) :: n_bodies ! Total number of bodies
integer(4) :: imping_body_grav ! Flag to consider (1) or not (0) gravity
                               ! components perpendicular &
                               ! to the normal vectors of the impingement
                               ! surfaces
integer(4) :: NumBEdges ! Number of convex edges
double precision :: eta ! 0.001 * Domain%h
double precision :: eta2 ! 0.01 * Domain%h * Domain%h
double precision :: doubleh ! 2.*Domain%h
double precision :: squareh ! Domain%h*Domain%h
double precision :: doublesquareh ! doubleh * doubleh
double precision :: cubich ! Domain%h*Domain%h*Domain%h
double precision :: unosuh ! 1./Domain%h
double precision :: unosusquareh ! 1./(Domain%h*Domain%h)
double precision :: dx_dxbodies ! Ratio between fluid particle and body
                                ! particle size 
! Indices of cells that must be considered around the current one 
! in subrutine "CalcVarLength"
integer(4),dimension(14,3) :: indicecelle 
integer(4),dimension(0:3,2) :: icoordp ! Pointer for coordinate location 2D or 
                                       ! 3D; = (/0,1,3,0,0,1,2,3/) initialized
                                       ! in the main program for compatibility
                                       ! with xlf90
character(len=8),parameter :: acode = "SPHERA  "
character(len=8),parameter :: version = "7.2   "
character(255) :: nomecaso, nomecas2
character(1),dimension(0:3) :: xyzlabel = (/ "T", "X", "Y", "Z" /)  
character(4),dimension(3) :: ncordlabel = (/ "    ", "(2D)", "(3D)" /)  
character(255) :: LicenseFile ! File name for license 
character(255) :: nomefilekill ! Killer file name 
character(10) :: exetype ! Type of machine: "windows" or "linux"
character(255) :: nomefileerr ! File name for error file in erosion model
character(8) :: modelloerosione ! type of erosion criterion (shields, mohr)
! "original" or "euristic"
character(len=8) :: dt_opt = "original" 
character(80),dimension(MAXTIT) :: title  
end module

