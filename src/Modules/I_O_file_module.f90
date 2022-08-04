!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: I_O_file_module            
! Description: Module for I/O units: 11-max_file_unit_booked
!-------------------------------------------------------------------------------
module I_O_file_module
! Error unit
integer(4) :: uerr = 0
! Input
integer(4) :: ninp = 11
! Log unit
integer(4) :: ulog = 12
#if (defined SOLID_BODIES) && (defined SPACE_3D)
! ".vtu" file
integer(4) :: uvtu = 13
#endif
! Time-dependent restart file to write
integer(4) :: nres = 21
! Time-dependent restart file to read
integer(4) :: nsav = 22
! Free surface
integer(4) :: nplb = 23
! Fluid front 
integer(4) :: nfro = 24
! Monitoring lines/points 
integer(4) :: ncpt = 25
! Lower fluid top height
integer(4) :: uzlft = 26
#ifdef SPACE_3D
! 2D output fields
integer(4) :: u2DO = 27
! Monitoring sections (output)
integer(4) :: uQse = 28
#endif
! Paraview
integer(4) :: unitvtk = 29
! Post-processing of interfaces
integer(4) :: uint = 30
! Input external file 
integer(4) :: ninp2 = 31
! Dummy file
integer(4) :: ndum = 32
#ifdef SOLID_BODIES
! Body dynamics time series
integer(4) :: ubod = 33
#endif
#ifdef SPACE_3D
! 2D output file of CLC class and z0
integer(4) :: uCLC = 34
! Time-independent restart file for z0 to read
integer(4) :: urz0 = 35
! Time-independent restart file for z0 to write
integer(4) :: usz0 = 36
#endif
! Killer file
integer(4) :: unitkill = 51
! Elapsed time
integer(4) :: unit_time_elapsed = 52 
! Error file for erosion model
integer(4) :: ueroerr = 55
! Surface mesh list for DB-SPH
integer(4) :: unit_file_list = 56
! Surface mesh files for DB-SPH
integer(4) :: unit_DBSPH_mesh = 57
! DB-SPH post-processing (selection of surface element IDs) to write the 
! surface element values     
integer(4) :: unit_dbsph_se_ID = 58
! DB-SPH post-processing (selection of a domain region) to write Force along 
! x-axis
integer(4) :: unit_dbsph_Fx = 59
! DB-SPH post-processing (selection of a domain region) to write the surface 
! element values
integer(4) :: unit_dbsph_se_reg = 60
#ifdef SPACE_3D
! Substation post-processing
integer(4) :: unit_substations = 61
#endif
! Do not use in this module unit file numbers larger than max_file_unit_booked: 
! they might be assigned locally to temporary files.
character(255), dimension(0:7) :: nomefile
end module
