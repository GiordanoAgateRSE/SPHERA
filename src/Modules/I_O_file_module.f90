!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: I_O_file_module            
! Description: Module for I/O.                     
!----------------------------------------------------------------------------------------------------------------------------------

! I/O units for files
module I_O_file_module
integer(4) :: nscr = 0 ! Screen 
integer(4) :: ninp = 11 ! Input
integer(4) :: nout = 12 ! Log 
integer(4) :: nres = 21 ! Results  
integer(4) :: nsav = 22 ! Restart 
integer(4) :: nplb = 23 ! Free surface 
integer(4) :: nfro = 24 ! Fluid front 
integer(4) :: ncpt = 25 ! Monitoring lines/points 
integer(4) :: unitvtk = 29 ! Paraview 
integer(4) :: ninp2 = 31 ! Input external file
integer(4) :: ndum = 32 ! Dummy file
integer(4) :: unitkill = 51 ! Killer file
integer(4) :: unit_time_elapsed = 52 ! Elapsed time
integer(4) :: uniterr = 55 ! Error file for erosion model
integer(4) :: unit_file_list = 56 ! Surface mesh list for DB-SPH 
integer(4) :: unit_DBSPH_mesh = 57 ! Surface mesh files for DB-SPH
integer(4) :: unit_dbsph_se_ID = 58 ! DB-SPH post-processing 
                                    ! (selection of surface element IDs) to 
                                    ! write the surface element values     
integer(4) :: unit_dbsph_Fx = 59 ! DB-SPH post-processing
                                 ! (selection of a domain region) to write 
                                 ! Force along x-axis
integer(4) :: unit_dbsph_se_reg = 60 ! DB-SPH post-processing
                                     ! (selection of a domain region) to 
                                     ! write the surface element values
character(255), dimension(0:7) :: nomefile
end module
