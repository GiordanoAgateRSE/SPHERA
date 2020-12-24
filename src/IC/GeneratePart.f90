!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: GeneratePart
! Description: Particle positions (initial conditions).
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine GeneratePart(IC_loop)
#elif defined SPACE_2D
subroutine GeneratePart
#endif
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
#ifdef SPACE_3D
use I_O_file_module
use I_O_diagnostic_module
#endif
!------------------------
! Declarations
!------------------------
implicit none
#ifdef SPACE_3D
integer(4),intent(in) :: IC_loop
#endif
integer(4) :: Nz,Mate,IsopraS,NumParticles,i,NumPartPrima,test_z
#ifdef SPACE_3D
integer(4) :: aux_factor,i_vertex,j_vertex,test_xy,test_xy_2,test_face,test_dam
integer(4) :: n_levels,nag_aux,alloc_stat,i_face,j,k,npi
double precision :: z_min,distance_hor,aux_scal,rnd
double precision :: aux_vec(3)
#endif
integer(4),dimension(SPACEDIM) :: Npps
#ifdef SPACE_3D
double precision,dimension(NPartZone) :: h_reservoir
#endif
double precision,dimension(SPACEDIM) :: MinOfMin,XminReset
#ifdef SPACE_3D
double precision,dimension(:),allocatable :: z_aux
#endif
double precision,dimension(:,:),allocatable :: Xmin,Xmax
#ifdef SPACE_3D
character(len=lencard)  :: nomsub = "GeneratePart"
type(TyParticle),dimension(:),allocatable :: pg_aux
#endif
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
      double precision :: dis1,dis2
      double precision :: normal(2)
   end subroutine point_inout_convex_non_degenerate_polygon
end interface
!------------------------
! Allocations
!------------------------
allocate(Xmin(SPACEDIM,NPartZone),Xmax(SPACEDIM,NPartZone))
!------------------------
! Initializations
!------------------------
call random_seed()
NumParticles = 0
test_z = 0
#ifdef SPACE_3D
if (IC_loop==2) then
   do Nz=1,NPartZone
      if (Partz(Nz)%IC_source_type==2) then
         test_z = 1
      endif
   enddo
endif
#endif
if (test_z==0) nag = 0
#ifdef SPACE_3D
do Nz=1,NPartZone
   if (Domain%tipo=="bsph") then
! In case of DB-SPH boundary treatment, there is a fictitious fluid reservoir
! top, which completes the kernel suppport (only for pre-processing)
      h_reservoir(Nz) = Partz(Nz)%H_res + 3.d0 * Domain%h
      else
         h_reservoir(Nz) = Partz(Nz)%H_res
   endif
enddo
#endif
MinOfMin = max_positive_number
!------------------------
! Statements
!------------------------
! Loop over the zones in order to set the initial minimum and maximum 
! coordinates of each zone and of the numerical domain
first_cycle: do Nz=1,NPartZone
! To search for the maximum and minimum "partzone" coordinates
   Partz(Nz)%coordMM(1:SPACEDIM,1) = max_positive_number
   Partz(Nz)%coordMM(1:SPACEDIM,2) = max_negative_number
   Xmin(1:SPACEDIM,Nz) = max_positive_number
   Xmax(1:SPACEDIM,Nz) = max_negative_number
! To search for the minimum and maximum coordinates of the current zone 
#ifdef SPACE_3D
      call FindFrame(Xmin,Xmax,Nz)
#elif defined SPACE_2D
         call FindLine(Xmin,Xmax,Nz)
#endif
! To evaluate the minimum and maximum coordinates of the zone
   do i=1,SPACEDIM
      if (Xmin(i,Nz)<Partz(Nz)%coordMM(i,1)) Partz(Nz)%coordMM(i,1) = Xmin(i,Nz)
      if (Xmax(i,Nz)>Partz(Nz)%coordMM(i,2)) Partz(Nz)%coordMM(i,2) = Xmax(i,Nz)
! To evaluate the overall minimum among the zones
      MinOfMin(i) = min(MinOfMin(i),Xmin(i,Nz))
   enddo
enddo first_cycle
! Loop over the zone in order to set the particle locations
second_cycle: do Nz=1,NPartZone
! To skip the boundaries not having types equal to "perimeter" or "pool"
   if ((Partz(Nz)%tipo/="peri").and.(Partz(Nz)%tipo/="pool")) cycle second_cycle
! "perimeter" or "pool" conditions are set
   Mate = Partz(Nz)%Medium
   Partz(Nz)%limit(1) = NumParticles + 1
   if ((Partz(Nz)%tipo=="peri").and.(Partz(Nz)%IC_source_type==2)) then
#ifdef SPACE_3D
! During the first IC loop
      if (IC_loop==1) then
! Update number of points/vertices of the zone
         Partz(Nz)%npoints = Partz(Nz)%ID_last_vertex -                        &
                             Partz(Nz)%ID_first_vertex + 1
         aux_factor = nint(Partz(Nz)%dx_CartTopog / Domain%dx)
! Allocation of the auxiliary array z_aux
         if (.not.allocated(z_aux)) then
! The number of points is increased due to the eventual presence of the 
! fictitious vertex n.1
            allocate(z_aux(Partz(Nz)%npoints+1),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*) 'Allocation of the auxiliary array z_aux ',       &
                 'failed; the program stops here (subroutine GeneratePart). '
               stop
            endif
         endif
! Loops over Cartesian topography points
         z_min = max_positive_number
         do i_vertex=Partz(Nz)%ID_first_vertex,Partz(Nz)%ID_last_vertex
            z_min = min(z_min,(Vertice(3,i_vertex)))
         enddo
         nag_aux = Domain%nag_aux
! Allocation of the auxiliary array pg_aux
         if (.not.allocated(pg_aux)) then
            allocate(pg_aux(nag_aux),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*) 'Allocation of the auxiliary array pg_aux',       &
                 'failed; the program stops here (subroutine GeneratePart). '
               stop
            endif
         endif
         NumParticles = NumParticles + 1
! Loops over Cartesian topography points
         do i_vertex=Partz(Nz)%ID_first_vertex,Partz(Nz)%ID_last_vertex
! Check if the vertex is inside the plan_reservoir   
            call point_inout_convex_non_degenerate_polygon(                    &
               Vertice(1:2,i_vertex),Partz(Nz)%plan_reservoir_points,          &
               Partz(Nz)%plan_reservoir_pos(1,1:2),                            &
               Partz(Nz)%plan_reservoir_pos(2,1:2),                            &
               Partz(Nz)%plan_reservoir_pos(3,1:2),                            &
               Partz(Nz)%plan_reservoir_pos(4,1:2),                            &
               Partz(Nz)%plan_reservoir_pos(4,1:2),                            &
               Partz(Nz)%plan_reservoir_pos(4,1:2),test_xy)
            if (test_xy==1) then
! Set the minimum topography height among the closest 9 points
               z_aux(i_vertex) = max_positive_number
! Suggestion.
! This loop could be optimized not to consider all the DTM points.
               do j_vertex=1,Partz(Nz)%npoints
                  distance_hor = dsqrt((Vertice(1,i_vertex) -                  &
                                 Vertice(1,j_vertex)) ** 2 +                   &
                                 (Vertice(2,i_vertex) -                        &
                                 Vertice(2,j_vertex)) ** 2)
                  if (distance_hor<=(1.5*Partz(Nz)%dx_CartTopog))              &
                     z_aux(i_vertex) = min(z_aux(i_vertex),                    &
                                           (Vertice(3,j_vertex)))
               enddo
! Generating daughter points of the Cartesian topography, according to dx
! (maybe the function "ceiling" would be better than "int" here)
               n_levels = int((h_reservoir(Nz) - z_aux(i_vertex)) / Domain%dx)
               do i=1,aux_factor
                  do j=1,aux_factor
                     do k=1,n_levels
! Setting coordinates                       
! We avoid particles to be located on the very diagonals of topography 
! (otherwise some particles could be not removed even if below the topography; 
! it is unclear why this rare eventuality would happen, 
! but no problem remains anymore)
                        pg_aux(NumParticles)%coord(1) = Vertice(1,i_vertex)    &
                                                        - aux_factor / 2.d0    &
                                                        * Domain%dx + (i -     &
                                                        1 + 0.501d0) *         &
                                                        Domain%dx
                        pg_aux(NumParticles)%coord(2) = Vertice(2,i_vertex)    &
                                                        - aux_factor / 2.d0    &
                                                        * Domain%dx + (j -     &
                                                        1 + 0.5d0) *           &
                                                        Domain%dx
                        pg_aux(NumParticles)%coord(3) = (h_reservoir(Nz) -     &
                                                        Domain%dx / 2.d0) -    &
                                                        (k - 1) * Domain%dx
! Test if the particle is below the reservoir
! Loop over boundaries
                        test_z = 0
                        test_face = 0
!$omp parallel do default(none)                                                &
!$omp shared(NumFacce,Tratto,BoundaryFace,Partz,Nz,pg_aux,test_face,test_z)    &
!$omp shared(NumParticles)                                                     &
!$omp private(i_face,aux_vec,aux_scal,test_xy)                         
                        do i_face=1,NumFacce
                           if (test_face==1) cycle
                           if (Tratto(BoundaryFace(i_face)%stretch)%zone==     &
                              Partz(Nz)%Car_top_zone) then  
! Test if the point lies inside the plan projection of the face     
                              call point_inout_convex_non_degenerate_polygon(  &
                                 pg_aux(NumParticles)%coord(1:2),              &
                                 BoundaryFace(i_face)%nodes,                   &
                                 BoundaryFace(i_face)%Node(1)%GX(1:2),         &
                                 BoundaryFace(i_face)%Node(2)%GX(1:2),         &
                                 BoundaryFace(i_face)%Node(3)%GX(1:2),         &
                                 BoundaryFace(i_face)%Node(4)%GX(1:2),         &
                                 BoundaryFace(i_face)%Node(4)%GX(1:2),         &
                                 BoundaryFace(i_face)%Node(4)%GX(1:2),         &                                 
                                 test_xy)
                              if (test_xy==1) then
! No need for a critical section to update test_face: only one face/process 
! is interested (no conflict); in particular cases there could be a conflict, 
! but no face has a preference and test_face cannot come back to zero 
! (no actual problem)
                                 test_face = 1
! Test if the point is invisible (test_z=1, below the topography) or visible 
! (test_z=0) to the face
                                 aux_vec(:) = pg_aux(NumParticles)%coord(:)    &
                                    - BoundaryFace(i_face)%Node(1)%GX(:)
                                 aux_scal = dot_product                        &
                                    (BoundaryFace(i_face)%T(:,3),aux_vec)
! No need for a critical section to update test_z: in particular cases there 
! could be a conflict, but the "if" condition would always provide the same 
! result
                                 if (aux_scal<=0.d0) test_z=1
                              endif
                           endif
                        enddo
!$omp end parallel do
! Check the presence of a dam zone
                        test_dam = 0
                        if (Partz(Nz)%dam_zone_ID>0) then
! Test if the point lies inside the plan projection of the dam zone
                           call point_inout_convex_non_degenerate_polygon(     &
                              pg_aux(NumParticles)%coord(1:2),                 &
                              Partz(Nz)%dam_zone_n_vertices,                   &
                              Partz(Nz)%dam_zone_vertices(1,1:2),              &
                              Partz(Nz)%dam_zone_vertices(2,1:2),              &
                              Partz(Nz)%dam_zone_vertices(3,1:2),              &
                              Partz(Nz)%dam_zone_vertices(4,1:2),              &
                              Partz(Nz)%dam_zone_vertices(4,1:2),              &
                              Partz(Nz)%dam_zone_vertices(4,1:2),test_xy)
                           if (test_xy==1) then
                              test_face = 0
! Loop on the faces of the dam in the dam zone
!$omp parallel do default(none)                                                &
!$omp shared(NumFacce,test_dam,Tratto,BoundaryFace,Partz,Nz,pg_aux)            &
!$omp shared(NumParticles,test_face)                                           &
!$omp private(i_face,aux_vec,aux_scal,test_xy_2)
                              do i_face=1,NumFacce
                                 if (test_face==1) cycle
! Check if the face belongs to the dam_zone boundary and in particular to its 
! top
if ((Tratto(BoundaryFace(i_face)%stretch)%zone==Partz(Nz)%dam_zone_ID).and.    &
                                       (BoundaryFace(i_face)%T(3,3)<0.)) then 
! Test if the particle horizontal coordinates lie inside the horizontal 
! projection of the face
 call point_inout_convex_non_degenerate_polygon(                               &
                                       pg_aux(NumParticles)%coord(1:2),        &
                                       BoundaryFace(i_face)%nodes,             &
                                       BoundaryFace(i_face)%Node(1)%GX(1:2),   &
                                       BoundaryFace(i_face)%Node(2)%GX(1:2),   &
                                       BoundaryFace(i_face)%Node(3)%GX(1:2),   &
                                       BoundaryFace(i_face)%Node(4)%GX(1:2),   &
                                       BoundaryFace(i_face)%Node(4)%GX(1:2),   &
                                       BoundaryFace(i_face)%Node(4)%GX(1:2),   &
                                       test_xy_2)
                                    if (test_xy_2==1) then
! Test if the particle position is invisible (test_dam=1) or visible 
! (test_dam=0) to the dam. Visibility to the dam means invisibility to the 
! current dam top face. 
                                       aux_vec(:) =                           &
                                          pg_aux(NumParticles)%coord(:) -     &
                                          BoundaryFace(i_face)%Node(1)%GX(:)
                                       aux_scal = dot_product(                &
                                          BoundaryFace(i_face)%T(:,3),aux_vec)
! Note: even if a particle simultanously belongs to 2 faces test_dam will 
! provide the same results: no matter the face order
!$omp critical (GeneratePart_cs)
                                       test_face = 1
                                       if (aux_scal<0.d0) then 
                                          test_dam=0
                                          else
                                             test_dam = 1
                                       endif
!$omp end critical (GeneratePart_cs)
                                    endif   
                                 endif
                              enddo
!$omp end parallel do                      
                           endif   
                        endif
! Update fluid particle counter
! Check if the fluid particle is visible both to the bottom reservoir and the 
! eventual dam
                        if ((test_z==0).and.(test_dam==0)) then
                           NumParticles = NumParticles + 1 
! Check the storage for the reached number of fluid particles
                           if (NumParticles>nag_aux) then
                              call diagnostic(arg1=10,arg2=4,arg3=nomsub)
                           endif
                        endif
                     enddo
                  enddo
               enddo
            endif
         enddo
         Partz(Nz)%limit(2) = NumParticles - 1
         else
! Second IC loop
! Loops over fluid particles 
!$omp parallel do default(none)                                                &
!$omp shared(Domain,Nz,Mate,pg,Partz)                                          &
!$omp private(npi,rnd)
            do npi=Partz(Nz)%limit(1),Partz(Nz)%limit(2)
! To set particle parameters
               call SetParticleParameters(npi,Nz,Mate)
! To impose a white noise to particle positions
               if (Domain%RandomPos=='r') then
                  call random_number(rnd)
                  pg(npi)%coord(1) = pg(npi)%coord(1) + (2.d0 * rnd - 1.d0) *  &
                                     0.1d0 * Domain%dx
                  call random_number(rnd)
                  pg(npi)%coord(2) = pg(npi)%coord(2) + (2.d0 * rnd - 1.d0) *  &
                                     0.1d0 * Domain%dx
                  call random_number(rnd)
                  pg(npi)%coord(3) = pg(npi)%coord(3) + (2.d0 * rnd - 1.d0) *  &
                                     0.1d0 * Domain%dx
               endif
               pg(npi)%sect_old_pos(:) = pg(npi)%coord(:)
            enddo
!$omp end parallel do
            NumParticles = Partz(Nz)%limit(2)
      endif
#endif
      else
! To evaluate the number of particles Npps for the zone along all the 
! directions
         if ((Xmin(1,Nz)==max_positive_number).or.                             &
             (Xmax(1,Nz)==max_negative_number)) then
! The zone is declared but is not used
            Npps = -1
            else
               do i=1,SPACEDIM
                  NumPartPrima = nint((Xmin(i,Nz) - MinOfMin(i)) / Domain%dx)
                  Xminreset(i) = MinOfMIn(i) + NumPartPrima * Domain%dx
                  Npps(i) = nint((Xmax(i,Nz) - XminReset(i)) / Domain%dx)
               enddo
         endif
#ifdef SPACE_2D
! To force almost one particle if the domain width in a direction is smaller  
! than dx
            where (Npps==0) Npps = 1
#endif
! To consider the "pool" condition: set "Isopra=1" if it crosses the 
! "virtual" side
         if (Partz(Nz)%tipo=="pool") then
            Xmax(Partz(Nz)%ipool,Nz) = Partz(Nz)%pool
            IsopraS = 1
            else
               IsopraS = 0
         endif
! To set the particles in the zone
#ifdef SPACE_3D
            call SetParticles(Nz,Mate,XminReset,Npps,NumParticles,IsopraS)
#elif defined SPACE_2D
               call SetParticles(Nz,Mate,XminReset,Npps,NumParticles)
#endif
! To set the upper pointer for the particles in the nz-th zone
         Partz(Nz)%limit(2) = NumParticles
   endif
enddo second_cycle
#ifdef SPACE_3D
if (allocated(pg_aux)) then
! Allocation of the fluid particle array in the presence of at least one  
! reservoir extruded from topography       
   NumParticles = NumParticles - 1
   PARTICLEBUFFER = int(NumParticles * Domain%COEFNMAXPARTI) + 1
   allocate(pg(PARTICLEBUFFER),stat=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,'(1x,a,i2)')                                                 &
         "    Array PG not allocated. Error code: ",alloc_stat
      call diagnostic(arg1=4,arg3=nomsub)
      else
         write(ulog,'(1x,a)') "    Array PG successfully allocated "
         pg(:) = PgZero
   endif
! Loop over the auxiliary particle array
!$omp parallel do default(none) shared(pg,pg_aux,NumParticles) private(npi)
   do npi=1,NumParticles
! Copy the fluid particle array from the corresponding auxiliary array
      pg(npi) = pg_aux(npi)
   enddo
!$omp end parallel do
! Deallocation of the auxiliary array z_aux
   deallocate(z_aux,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Deallocation of the auxiliary array z_aux ',              &
         'failed; the program stops here (subroutine GeneratePart). '
      stop
   endif
! Deallocation of the auxiliary array pg_aux
   deallocate(pg_aux,STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(ulog,*) 'Deallocation of the auxiliary array pg_aux ',             &
         'failed; the program stops here (subroutine GeneratePart). '
      stop
   endif
endif
#endif
test_z = 0
#ifdef SPACE_3D
if (IC_loop==2) then
   do Nz=1,NPartZone
      if (Partz(Nz)%IC_source_type==2) test_z = 1
   enddo
endif
#endif
if (test_z==0) nag = NumParticles
nagpg = nag
!------------------------
! Deallocations
!------------------------
deallocate(Xmin,Xmax)
return
end subroutine GeneratePart
