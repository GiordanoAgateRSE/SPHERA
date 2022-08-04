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
! Program unit: electrical_substations    
! Description: output assessment and writing for electrical substations (only 
!              in 3D). Output (time series): Probability of an Outage Start 
!              (POS), Expected Outage Status (EOS), Expected Outage Time (EOT, 
!              time series update of its expected scalar value), Damage to the 
!              Substation Beyond Outage (Dsub, time series update of its 
!              expected scalar value), Substation Vulnerability (Vul). 
!              Output depends on Ysub (spatial average of the fluid/mixture 
!              depth at the DEM grid points, within the substation polygon). 
!              Ysub alternatively depends on two estimations of the water 
!              depth: the first (i_fil==1) is a raw estimation; the latter 
!              (i_fil==2) filters the atomization and the wave breaking effects.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine electrical_substations
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
integer(4) :: test,i_zone,i_vertex,i_sub,i_aux,GridColumn,aux_integer,DEM_points
integer(4) :: alloc_stat,i_fil
double precision :: Y_step,Y_filt_step
character(255) :: nomefile_substations
double precision :: aux_vec(3)
! Expected Outage Status
logical,dimension(:,:),allocatable :: EOS
! Probability of an outage start
double precision,dimension(:,:),allocatable :: POSsub
! Substation fluid/mixture depth
double precision,dimension(:,:),allocatable :: Ysub
! Damage to the substation
double precision,dimension(:,:),allocatable :: Dsub
! Substation vulnerability
double precision,dimension(:,:),allocatable :: Vul
double precision :: pos(6,3)
integer(4),dimension(:,:),allocatable :: aux_array
integer(4),external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine point_inout_convex_non_degenerate_polygon(point,n_sides,         &
                                                        point_pol_1,           &
                                                        point_pol_2,           &
                                                        point_pol_3,           &
                                                        point_pol_4,           &
                                                        point_pol_5,           &
                                                        point_pol_6,test)
      implicit none
      integer(4),intent(in) :: n_sides
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_convex_non_degenerate_polygon
   subroutine point_inout_hexagon(point,point_pol_1,point_pol_2,               &
                                  point_pol_3,point_pol_4,point_pol_5,         &
                                  point_pol_6,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_hexagon
   subroutine point_inout_pentagon(point,point_pol_1,point_pol_2,              &
                                   point_pol_3,point_pol_4,point_pol_5,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_pentagon
   subroutine point_inout_quadrilateral(point,point_pol_1,point_pol_2,         &
                                        point_pol_3,point_pol_4,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_quadrilateral
   subroutine area_hexagon(P1,P2,P3,P4,P5,P6,area)
      implicit none
      double precision,intent(IN) :: P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
      double precision,intent(inout) :: area
   end subroutine area_hexagon  
   subroutine area_pentagon(P1,P2,P3,P4,P5,area)
      implicit none
      double precision,intent(in) :: P1(3),P2(3),P3(3),P4(3),P5(3)
      double precision,intent(inout) :: area
   end subroutine area_pentagon
   subroutine area_quadrilateral(P1,P2,P3,P4,area)
      implicit none
      double precision,intent(in) :: P1(3),P2(3),P3(3),P4(3)
      double precision,intent(inout) :: area
   end subroutine area_quadrilateral
   subroutine area_triangle(P1,P2,P3,area,normal)
      implicit none
      double precision,intent(in) :: P1(3),P2(3),P3(3)
      double precision,intent(out) :: area
      double precision,intent(out) :: normal(3)
   end subroutine area_triangle
end interface
!------------------------
! Allocations
!------------------------
if (.not.allocated(EOS)) then
   allocate(EOS(substations%n_sub,2),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "EOS" failed in the subroutine ',           &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (.not.allocated(POSsub)) then
   allocate(POSsub(substations%n_sub,2),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "POSsub" failed in the subroutine ',        &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (.not.allocated(Ysub)) then
   allocate(Ysub(substations%n_sub,2),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "Ysub" failed in the subroutine ',          &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (.not.allocated(Dsub)) then
   allocate(Dsub(substations%n_sub,2),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "Dsub" failed in the subroutine ',          &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (.not.allocated(Vul)) then
   allocate(Vul(substations%n_sub,2),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Allocation of "Vul" failed in the subroutine ',           &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
!------------------------
! Initializations
!------------------------
write(nomefile_substations,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),     &
   '_substations_',on_going_time_step,".txt"
EOS(:,:) = .false.
POSsub(:,:) = 0.d0
Ysub(:,:) = 0.d0
Dsub(:,:) = 0.d0
Vul(:,:) = 0.d0
!------------------------
! Statements
!------------------------
open(unit_substations,file=nomefile_substations,status="unknown",              &
     form="formatted")
if (.not.allocated(substations%sub(1)%DEMvert)) then
! First time step in any case (both with/without restart)
! .txt file creation and headings
   write(unit_substations,*) "Electrical substations "
   write(unit_substations,                                                     &
'((7x,a),(1x,a),(4x,a),(8x,a),(4x,a),(10x,a),(7x,a),2(11x,a),(8x,a),(11x,a),(3x,a))'&
      ) " Time(s)"," ID_substation"," Val(euros)"," ID_fil"," area(m**2)",     &
      " Y(m)"," Ymax(m)"," POS"," EOS"," EOT(s)"," Vul"," Dsub(euros)"
   flush(unit_substations)
! Association of the DEM points with the substations: start.
! Loop over the zones
   do i_zone=1,NPartZone
      if (Partz(i_zone)%ID_first_vertex_sel>0) then
         DEM_points = Partz(i_zone)%ID_last_vertex_sel -                       &
                      Partz(i_zone)%ID_first_vertex_sel + 1
! Allocation and initialization of the auxiliary array for the DEM points 
! associated with the substations.
         if (.not.allocated(aux_array)) then
            allocate(aux_array(substations%n_sub,DEM_points),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(uerr,*) 'Allocation of "aux_array" failed in the ',       &
                  'subroutine "electrical_substations". The execution ',       &
                  'terminates here.'
               stop
               else
                  write(ulog,*) 'Allocation of "aux_array" in the subroutine', &
                     '"electrical_substations" is successfully completed.'
            endif
            aux_array(:,:) = 0
         endif
! Loop over the DEM/DTM vertices
!$omp parallel do default(none)                                                &
!$omp shared(Partz,i_zone,Vertice,substations,uerr,aux_array)                  &
!$omp private(i_vertex,pos,i_sub,test)
         do i_vertex=Partz(i_zone)%ID_first_vertex_sel,                        &
            Partz(i_zone)%ID_last_vertex_sel
            pos(1,1) = Vertice(1,i_vertex)
            pos(1,2) = Vertice(2,i_vertex)
! Loop over the substations
            do i_sub=1,substations%n_sub
              select case (substations%sub(i_sub)%n_vertices)
                  case(3)
                     call point_inout_convex_non_degenerate_polygon(pos(1,1:2),&
                             3,substations%sub(i_sub)%vert(1,1:2),             &
                             substations%sub(i_sub)%vert(2,1:2),               &
                             substations%sub(i_sub)%vert(3,1:2),               &
                             substations%sub(i_sub)%vert(4,1:2),               &
                             substations%sub(i_sub)%vert(4,1:2),               &
                             substations%sub(i_sub)%vert(4,1:2),test)
                  case(4)
                     call point_inout_quadrilateral(pos(1,1:2),                &
                             substations%sub(i_sub)%vert(1,1:2),               &
                             substations%sub(i_sub)%vert(2,1:2),               &
                             substations%sub(i_sub)%vert(3,1:2),               &
                             substations%sub(i_sub)%vert(4,1:2),test)
                  case(5)
                     call point_inout_pentagon(pos(1,1:2),                     &
                             substations%sub(i_sub)%vert(1,1:2),               &
                             substations%sub(i_sub)%vert(2,1:2),               &
                             substations%sub(i_sub)%vert(3,1:2),               &
                             substations%sub(i_sub)%vert(4,1:2),               &
                             substations%sub(i_sub)%vert(5,1:2),test)
                  case(6)
                     call point_inout_hexagon(pos(1,1:2),                      &
                             substations%sub(i_sub)%vert(1,1:2),               &
                             substations%sub(i_sub)%vert(2,1:2),               &
                             substations%sub(i_sub)%vert(3,1:2),               &
                             substations%sub(i_sub)%vert(4,1:2),               &
                             substations%sub(i_sub)%vert(5,1:2),               &
                             substations%sub(i_sub)%vert(6,1:2),test)
                  case default
                     write(uerr,*) 'Error in defining the polygons of the ',   &
                        'electrical substations (subroutine ',                 &
                        '"electrical_substations"). The execution terminates ',&
                        'here. '
                     stop
               endselect
               if (test==1) then
                  substations%sub(i_sub)%n_DEM_vertices =                      &
                     substations%sub(i_sub)%n_DEM_vertices + 1
                  aux_array(i_sub,substations%sub(i_sub)%n_DEM_vertices) =     &
                     i_vertex
               endif
            enddo
         enddo
!$omp end parallel do
         exit
      endif
   enddo
! Association of the DEM points with the substations: end.
! Loop over the substations
!$omp parallel do default(none)                                                &
!$omp shared(substations,aux_array,uerr,ulog)                                  &
!$omp private(i_sub,aux_vec,alloc_stat,pos)
   do i_sub=1,substations%n_sub
! Allocation and assignment for the array of the substation DEM vertices
      if (.not.allocated(substations%sub(i_sub)%DEMvert)) then
allocate(substations%sub(i_sub)%DEMvert(substations%sub(i_sub)%n_DEM_vertices) &
            ,STAT=alloc_stat)
         if (alloc_stat/=0) then
            write(uerr,*) 'Allocation of "substations%sub(',i_sub,             &
               ')%DEMvert" in the subroutine "electrical_substations" failed', &
               '; the execution terminates here.'
            stop
            else
               write(ulog,*) 'Allocation of "substations%sub(',i_sub,          &
                  ')%DEMvert" in the subroutione "electrical_substations" is ',&
                  'successfully completed.'
         endif
      endif
      substations%sub(i_sub)%DEMvert(1:substations%sub(i_sub)%n_DEM_vertices) =&
         aux_array(i_sub,1:substations%sub(i_sub)%n_DEM_vertices)
! Substation area: start.
      pos(1:6,1:2) = substations%sub(i_sub)%vert(1:6,1:2)
      pos(1:6,3) = 0.d0
      select case (substations%sub(i_sub)%n_vertices)
         case(3)
            call area_triangle(pos(1,1:3),pos(2,1:3),pos(3,1:3),               &
                               substations%sub(i_sub)%area,aux_vec)
         case(4)
            call area_quadrilateral(pos(1,1:3),pos(2,1:3),pos(3,1:3),          &
                                    pos(4,1:3),substations%sub(i_sub)%area)
         case(5)
            call area_pentagon(pos(1,1:3),pos(2,1:3),pos(3,1:3),pos(4,1:3),    &
                               pos(5,1:3),substations%sub(i_sub)%area)
         case(6)
            call area_hexagon(pos(1,1:3),pos(2,1:3),pos(3,1:3),pos(4,1:3),     &
                              pos(5,1:3),pos(6,1:3),substations%sub(i_sub)%area)
         case default
                  write(uerr,*) 'Error in defining the areas of the ',         &
                     'electrical substations (subroutine ',                    &
                     '"electrical_substations"). The execution terminates here.'
                  stop
      endselect
! Substation area: end.
! Substation value
      select case (substations%sub(i_sub)%type_ID)
         case(1)
            substations%sub(i_sub)%Val = 41000.d0 
         case(2)
            substations%sub(i_sub)%Val = 16400.d0
         case(3)
            substations%sub(i_sub)%Val = 8200.d0
         case default
            write(uerr,*) 'Error in defining the values of the electrical ',   &
               'substations (subroutine "electrical_substations"). The ',      &
               'execution terminates here. '
               stop
      endselect
   enddo
!$omp end parallel do
   else
! Other time steps
! Loop over the substations
!$omp parallel do default(none)                                                &
!$omp shared(substations,Vertice,Grid,Z_fluid_step,EOS,Vul,Dsub,Ysub,POSsub)   &
!$omp private(i_sub,i_aux,pos,GridColumn,Y_step,Y_filt_step,i_fil)
      do i_sub=1,substations%n_sub
! Spatial average of the fluid/mixture depth at the DEM grid points, within the 
! substation polygon: start.
! Loop over the DEM points of the substation
         do i_aux=1,substations%sub(i_sub)%n_DEM_vertices
            pos(1,1) = Vertice(1,substations%sub(i_sub)%DEMvert(i_aux))
            pos(1,2) = Vertice(2,substations%sub(i_sub)%DEMvert(i_aux))
            pos(1,3) = Grid%extr(3,1) + 1.d-7
            GridColumn = ParticleCellNumber(pos(1,1:3))
            Y_step = max((Z_fluid_step(GridColumn,1) -                         &
                         Vertice(3,substations%sub(i_sub)%DEMvert(i_aux))),0.d0)
            Y_filt_step = max((Z_fluid_step(GridColumn,2) -                    &
               Vertice(3,substations%sub(i_sub)%DEMvert(i_aux))),0.d0)
            Ysub(i_sub,1) = Ysub(i_sub,1) + Y_step
            Ysub(i_sub,2) = Ysub(i_sub,2) + Y_filt_step
         enddo
         if (substations%sub(i_sub)%n_DEM_vertices>0) then
            Ysub(i_sub,:) = Ysub(i_sub,:) /                                    &
                            substations%sub(i_sub)%n_DEM_vertices
         endif
         do i_fil=1,2
! Spatial average of the fluid/mixture depth at the DEM grid points, within the 
! substation polygon: end.
! Update of the maximum value of the substation fluid/mixture depth
            substations%sub(i_sub)%Ymax(i_fil) = max(Ysub(i_sub,i_fil),        &
               substations%sub(i_sub)%Ymax(i_fil))
! Probability of an Outage Start: POS (time series update)
            if (Ysub(i_sub,i_fil)>=0.52d0) then
               POSsub(i_sub,i_fil) = 1.d0
               elseif (Ysub(i_sub,i_fil)>0.d0) then
                  POSsub(i_sub,i_fil) = 1.92d0 * Ysub(i_sub,i_fil)
                  else
                     POSsub(i_sub,i_fil) = 0.d0
            endif
! Expected Outage Status (EOS, time series update); Expected Outage Time (EOT, 
! scalar written as a time series update). Simplifying hypothesis: physical 
! simulated time shorter than 11h.
            substations%sub(i_sub)%POS_fsum(i_fil) =                           &
               substations%sub(i_sub)%POS_fsum(i_fil) +                        &
               max(POSsub(i_sub,i_fil) - 0.49d0,0.d0)
            if (substations%sub(i_sub)%POS_fsum(i_fil)>1.d-9) then
               EOS(i_sub,i_fil) = .true.
               substations%sub(i_sub)%EOT(i_fil) =                             &
                  substations%sub(i_sub)%EOT(i_fil) + substations%dt_out
            endif
! Vulnerability (scalar written as a time series update)
            if (substations%sub(i_sub)%Ymax(i_fil)>10.d0) then
               Vul(i_sub,i_fil) = 0.15d0
               elseif (substations%sub(i_sub)%Ymax(i_fil)>0.d0) then
                  Vul(i_sub,i_fil) = -1.22877d-6 *                             &
                                     substations%sub(i_sub)%Ymax(i_fil) ** 6 + &
                                     1.92478d-5 *                              &
                                     substations%sub(i_sub)%Ymax(i_fil) ** 5 + &
                                     8.4216d-6 *                               &
                                     substations%sub(i_sub)%Ymax(i_fil) ** 4 - &
                                     0.00119121d0 *                            &
                                     substations%sub(i_sub)%Ymax(i_fil) ** 3 + &
                                     0.00390726d0 *                            &
                                     substations%sub(i_sub)%Ymax(i_fil) ** 2 + &
                                     0.0170243d0 *                             &
                                     substations%sub(i_sub)%Ymax(i_fil)
                  else
                     Vul(i_sub,i_fil) = 0.d0
            endif
! Damage to the substation (scalar written as a time series update)
            Dsub(i_sub,i_fil) = Vul(i_sub,i_fil) * substations%sub(i_sub)%Val
         enddo
      enddo
!$omp end parallel do
! Output writing
! Loop over the substations
      do i_sub=1,substations%n_sub
         do i_fil=1,2
            aux_integer = 0
            if (EOS(i_sub,i_fil).eqv..true.) aux_integer = 1
            write(unit_substations,                                            &
               '(g15.5,i15,g15.5,i15,4(g15.5),i15,3(g15.5))') simulation_time, &
               i_sub,substations%sub(i_sub)%Val,i_fil,                         &
               substations%sub(i_sub)%area,Ysub(i_sub,i_fil),                  &
               substations%sub(i_sub)%Ymax(i_fil),POSsub(i_sub,i_fil),         &
               aux_integer,substations%sub(i_sub)%EOT(i_fil),Vul(i_sub,i_fil), &
               Dsub(i_sub,i_fil)
         enddo
      enddo
endif
close(unit_substations)
!------------------------
! Deallocations
!------------------------
if(allocated(aux_array)) then
   deallocate(aux_array,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "aux_array" in the subroutine ',          &
         '"electrical_substations" failed; the execution terminates here. '
      stop
   endif
endif
if (allocated(EOS)) then
   deallocate(EOS,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "EOS" failed in the subroutine ',         &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (allocated(POSsub)) then
   deallocate(POSsub,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "POSsub" failed in the subroutine ',      &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (allocated(Ysub)) then
   deallocate(Ysub,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "Ysub" failed in the subroutine ',        &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (allocated(Dsub)) then
   deallocate(Dsub,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "Dsub" failed in the subroutine ',        &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
if (allocated(Vul)) then
   deallocate(Vul,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(uerr,*) 'Deallocation of "Vul" failed in the subroutine ',         &
         '"electrical_substations". The execution terminates here.'
      stop
   endif
endif
return
end subroutine electrical_substations
#endif
