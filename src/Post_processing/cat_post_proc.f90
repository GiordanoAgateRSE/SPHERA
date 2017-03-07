!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: cat_post_proc                  
! Description: To concatenate the ".txt" output files and remove the original 
!              ones.              
!-------------------------------------------------------------------------------
subroutine cat_post_proc
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
call system("touch Body_dynamics.txt Body_dynamics_.tmp")
call system("cat Body_dynamics.txt *Body_dynamics_* > Body-dynamics.txt")
call system("rm -f *Body_dynamics*")
call system("mv Body-dynamics.txt Body_dynamics.txt")
call system("touch Body_particles.txt Body_particles_.tmp")
call system("cat Body_particles.txt *Body_particles_* > Body-particles.txt")
call system("rm -f *Body_particles*")
call system("mv Body-particles.txt Body_particles.txt")
call system("touch Q_sections.txt Q_sections_.tmp")
call system("cat Q_sections.txt *Q_sections_* > Q-sections.txt")
call system("rm -f *Q_sections*")
call system("mv Q-sections.txt Q_sections.txt")
call system("touch monitoring_lines.txt temp.cln")
call system("cat monitoring_lines.txt *.cln > monitoring-lines.txt")
call system("rm -f monitoring_lines.txt *.cln")
call system("mv monitoring-lines.txt monitoring_lines.txt")
call system("touch monitoring_points.txt temp.cpt")
call system("cat monitoring_points.txt *.cpt > monitoring-points.txt")
call system("rm -f monitoring_points.txt *.cpt")
call system("mv monitoring-points.txt monitoring_points.txt")
call system("touch blt_interfaces.txt blt_interfaces_.tmp")
call system("cat blt_interfaces.txt *blt_interfaces_* > blt-interfaces.txt")
call system("rm -f *blt_interfaces*")
call system("mv blt-interfaces.txt blt_interfaces.txt")
!------------------------
! Deallocations
!------------------------
return
end subroutine cat_post_proc

