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
! Program unit: sub_Q_sections                 
! Description: Writing flow rate at monitoring sections provided in input for the flow rate (only in 3D).           
!----------------------------------------------------------------------------------------------------------------------------------

subroutine sub_Q_sections
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,i_sect,test_intersection_point,test_inout,j
character(255) :: nomefile_Q_sections
double precision :: sign_dis_part_sect_old,sign_dis_part_sect
double precision :: P_intersection(3),aux_vec_1(3),aux_vec_2(3)
double precision :: P_intersection_loc(3)
integer(4), allocatable, dimension(:) :: n_particles 
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
allocate(n_particles(Q_sections%n_sect))
!------------------------
! Initializations
!------------------------
n_particles = 0
!------------------------
! Statements
!------------------------
! .txt file creation and headings (only at the beginning of the simulation)
write(nomefile_Q_sections,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),      &
   '_Q_sections_',it_corrente,".txt"
open(ncpt,file=nomefile_Q_sections,status="unknown",form="formatted")
if (it_corrente==1) then
! First step
   write (ncpt,*) "Flow rate values(m^3/s) "
   write (ncpt,'((7x,a),(4x,a),(5x,a),(6x,a),(5x,a),(3x,a),(5x,a),6(5x,a))')   &
      " Time(s)"," ID_section","fluid_type"," Q(m^3/s)"," Area(m^2)",          &
      " n_particles"," dt_out(s)"," norm_x(m)"," norm_y(m)"," norm_z(m)",      &
      " x1_sec(m)"," y1_sec(m)"," z1_sec(m)"
   flush(ncpt)
! Initializing fluid particle old coordinates 
!$omp parallel do default(none) shared(nag,Q_sections,pg) private(npi)
   do npi=1,nag
      pg(npi)%sect_old_pos(:) = pg(npi)%coord(:)
   enddo  
!$omp end parallel do
! Initializing section local axis and vertex local coordinates
!$omp parallel do default(none) shared(Q_sections) private(i_sect)
   do i_sect=1,Q_sections%n_sect
! Local axis 1: Vertex2-vertex1
      Q_sections%section(i_sect)%loc_axis(:,1) =                               &
         Q_sections%section(i_sect)%vertex(2,:) -                              &
         Q_sections%section(i_sect)%vertex(1,:)
      Q_sections%section(i_sect)%loc_axis(:,1) =                               &
         Q_sections%section(i_sect)%loc_axis(:,1) /                            &
         dsqrt(dot_product(Q_sections%section(i_sect)%loc_axis(:,1),           &
         Q_sections%section(i_sect)%loc_axis(:,1)))
! Local axis 2: opposite of the vector product of Local axis 1 and the normal 
      call Vector_Product(Q_sections%section(i_sect)%loc_axis(:,1),            &
         Q_sections%section(i_sect)%normal,                                    &
         Q_sections%section(i_sect)%loc_axis(:,2),3)                           
      Q_sections%section(i_sect)%loc_axis(:,2) = -                             &
         Q_sections%section(i_sect)%loc_axis(:,2)
! Local axis 3: normal vector
      Q_sections%section(i_sect)%loc_axis(:,3) =                               &
         Q_sections%section(i_sect)%normal(:)
      Q_sections%section(i_sect)%vertex_loc(1,:) = 0.0d0
      call reference_system_change(Q_sections%section(i_sect)%vertex(2,:),     &
         Q_sections%section(i_sect)%vertex(1,:),                               &
         Q_sections%section(i_sect)%loc_axis,                                  &
         Q_sections%section(i_sect)%vertex_loc(2,:)) 
      call reference_system_change(Q_sections%section(i_sect)%vertex(3,:),     &
         Q_sections%section(i_sect)%vertex(1,:),                               &
         Q_sections%section(i_sect)%loc_axis,                                  &
         Q_sections%section(i_sect)%vertex_loc(3,:)) 
      if (Q_sections%section(i_sect)%n_vertices==4) then
         call reference_system_change(Q_sections%section(i_sect)%vertex(4,:),  &
            Q_sections%section(i_sect)%vertex(1,:),                            &
            Q_sections%section(i_sect)%loc_axis,                               &
            Q_sections%section(i_sect)%vertex_loc(4,:))
      endif
! Allocating and initializing the flow rate vector
      allocate(Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1))
      Q_sections%section(i_sect)%flow_rate(:) = 0.0d0
   enddo
!$omp end parallel do 
   else
! Other steps     
! Loop over the fluid particles
!$omp parallel do default(none) shared(nag,Q_sections,pg,n_particles,Med)      &
!$omp private(npi,i_sect,sign_dis_part_sect_old,sign_dis_part_sect)            &
!$omp private(test_intersection_point,P_intersection,test_inout,aux_vec_1)     &
!$omp private(aux_vec_2,P_intersection_loc)
      do npi=1,nag
! Loop over the flow rate monitoring sections 
         do i_sect=1,Q_sections%n_sect
! To check if a fluid particle has crossed the plane of the section
            aux_vec_1 = pg(npi)%sect_old_pos -                                 &
                        Q_sections%section(i_sect)%vertex(1,:)
            aux_vec_2 = pg(npi)%coord - Q_sections%section(i_sect)%vertex(1,:)
            sign_dis_part_sect_old = dot_product                               &
               (aux_vec_1,Q_sections%section(i_sect)%normal)
            sign_dis_part_sect = dot_product(aux_vec_2,                        &
                                 Q_sections%section(i_sect)%normal)
            if ((sign_dis_part_sect*sign_dis_part_sect_old)<0.d0) then
! To check if the intersection between the line (defined by the old and new
! particle position) and the plane (containing the section) is a point            
               call line_plane_intersection(pg(npi)%sect_old_pos,pg(npi)%coord,&
                  Q_sections%section(i_sect)%vertex(1,:),                      &
                  Q_sections%section(i_sect)%vertex(2,:),                      &
                  Q_sections%section(i_sect)%vertex(3,:),                      &
                  test_intersection_point,P_intersection)
               if (test_intersection_point==1) then
! Local coordinates
                  call reference_system_change(P_intersection,                 &
                     Q_sections%section(i_sect)%vertex(1,:),                   &
                     Q_sections%section(i_sect)%loc_axis,P_intersection_loc) 
! Check if the intersection point lies within the section
                  call point_inout_polygone(P_intersection_loc,                &
                     Q_sections%section(i_sect)%n_vertices,                    &
                     Q_sections%section(i_sect)%vertex_loc(1,:),               &
                     Q_sections%section(i_sect)%vertex_loc(2,:),               &
                     Q_sections%section(i_sect)%vertex_loc(3,:),               &
                     Q_sections%section(i_sect)%vertex_loc(4,:),test_inout)
                  if (test_inout==1) then
!$omp critical (omp_Q_Sections)                 
! To update the number of crossed particles and the flowing mass                  
                     n_particles(i_sect) = n_particles(i_sect) + 1
                     Q_sections%section(i_sect)%flow_rate(pg(npi)%imed) =      &
                        Q_sections%section(i_sect)%flow_rate(pg(npi)%imed) +   &
                        pg(npi)%mass/Med(pg(npi)%imed)%den0                    &
                        * (sign_dis_part_sect)/abs(sign_dis_part_sect)
!$omp end critical (omp_Q_Sections)
                  endif
               endif
            endif
         enddo
! To update the "old coordinates"    
         pg(npi)%sect_old_pos(:) = pg(npi)%coord(:)
      enddo  
!$omp end parallel do 
! Loop over the flow rate monitoring sections 
      do i_sect=1,Q_sections%n_sect
         do j=1,Q_sections%n_fluid_types
! To compute the volume to flow rate per fluid type: flowing mass by time 
            Q_sections%section(i_sect)%flow_rate(j) =                          &
               Q_sections%section(i_sect)%flow_rate(j) / Q_sections%dt_out
! Writing the flow rate on a ".txt" file    
            write(ncpt,'((f14.6,1x),2(i14,1x),2(f14.6,1x),(i14,1x),7(f14.6,1x))'&
               )tempo,i_sect,j,Q_sections%section(i_sect)%flow_rate(j),        &
               Q_sections%section(i_sect)%area,n_particles(i_sect),            &
               Q_sections%dt_out,Q_sections%section(i_sect)%normal(:),         &
               Q_sections%section(i_sect)%vertex(1,:)
! Update of the global volume flow rate 
            Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1) = &
               Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1)&
               + Q_sections%section(i_sect)%flow_rate(j)       
! Zeroing the flow rate per fluid type    
            Q_sections%section(i_sect)%flow_rate(j) = 0.0d0
         enddo
! Writing the global flow rate on a ".txt" file  
         j = Q_sections%n_fluid_types + 1
         write(ncpt,'((f14.6,1x),2(i14,1x),2(f14.6,1x),(i14,1x),7(f14.6,1x))') &
            tempo,i_sect,j,Q_sections%section(i_sect)%flow_rate(j),            &
            Q_sections%section(i_sect)%area,n_particles(i_sect),               &
            Q_sections%dt_out,Q_sections%section(i_sect)%normal(:),            &
            Q_sections%section(i_sect)%vertex(1,:) 
! Zeroing the flow rate per fluid type    
         Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1) =    &
            0.0d0    
      enddo
endif 
close (ncpt)
deallocate(n_particles) 
!------------------------
! Deallocations
!------------------------
return
end subroutine sub_Q_sections

