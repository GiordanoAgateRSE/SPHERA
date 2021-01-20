!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
logical :: z_min_flag
integer(4),intent(in) :: IC_loop
#endif
integer(4) :: Nz,Mate,IsopraS,NumParticles,i,NumPartPrima,test_z
#ifdef SPACE_3D
integer(4) :: aux_factor,i_vertex,test_xy,test_dam,size_aux
integer(4) :: n_levels,nag_aux,alloc_stat,j,k,npi
double precision :: max_rnd_eps
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
interface
   subroutine z_min_max_DEM_DTM_9p_stencil(min_flag,i_zone,i_vertex,z_aux)
      implicit none
      logical,intent(in) :: min_flag
      integer(4),intent(in) :: i_zone,i_vertex
      double precision,intent(inout) :: z_aux
      integer(4) :: j_vertex
      double precision :: distance_hor
   end subroutine z_min_max_DEM_DTM_9p_stencil
end interface
interface
   subroutine particle_position_extrusion(i_vertex,aux_factor,ii,jj,kk,        &
      DEM_zone,zmax,test_z,pos)
      implicit none
      integer(4),intent(in) :: i_vertex,aux_factor,ii,jj,kk,DEM_zone
      double precision,intent(in) :: zmax
      integer(4),intent(inout) :: test_z
      double precision,intent(inout) :: pos(3)
      integer(4) :: test_face,i_face,test_xy
      double precision :: aux_scal
      double precision :: aux_vec(3)
   end subroutine particle_position_extrusion
end interface
interface
   subroutine particles_in_out_dams(fluid_zone,pos,test_dam)
      implicit none
      integer(4),intent(in) :: fluid_zone
      double precision, intent(in) :: pos(3)
      integer(4),intent(out) :: test_dam
      integer(4) :: test_xy,test_face,i_face,test_xy_2
      double precision :: aux_scal
      double precision :: aux_vec(3)
   end subroutine particles_in_out_dams
end interface
interface
   subroutine pos_plus_white_noise(max_rnd_eps,pos)
      implicit none
      double precision,intent(in) :: max_rnd_eps
      double precision,intent(inout) :: pos(3)
      double precision :: rnd(3)
   end subroutine pos_plus_white_noise
end interface
!------------------------
! Allocations
!------------------------
allocate(Xmin(SPACEDIM,NPartZone),Xmax(SPACEDIM,NPartZone))
!------------------------
! Initializations
!------------------------
NumParticles = 0
test_z = 0
#ifdef SPACE_3D
z_min_flag = .true.
if (IC_loop==2) then
   do Nz=1,NPartZone
      if (Partz(Nz)%IC_source_type==2) then
         test_z = 1
      endif
   enddo
endif
max_rnd_eps = 0.1d0
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
call domain_edges(MinOfMin,Xmin,Xmax)
! Loop over the zone in order to set the particle positions
do Nz=1,NPartZone
! To skip the "boundaries" different from "perimeter" or "pool"
   if ((Partz(Nz)%tipo/="peri").and.(Partz(Nz)%tipo/="pool")) cycle
! "perimeter" or "pool" conditions are set
   Mate = Partz(Nz)%Medium
   Partz(Nz)%limit(1) = NumParticles + 1
   if ((Partz(Nz)%tipo=="peri").and.(Partz(Nz)%IC_source_type==2)) then
#ifdef SPACE_3D
! During the first IC loop
      if (IC_loop==1) then
! Update number of points/vertices of the zone
         size_aux = Partz(Nz)%ID_last_vertex_sel -                             &
                    Partz(Nz)%ID_first_vertex_sel + 1
         aux_factor = nint(Partz(Nz)%dx_CartTopog / Domain%dx)
! Allocation of the auxiliary array z_aux
         if (.not.allocated(z_aux)) then
! The number of points is increased due to the eventual presence of the 
! fictitious vertex n.1
            allocate(z_aux(size_aux+1),STAT=alloc_stat)
            if (alloc_stat/=0) then
               write(ulog,*) 'Allocation of the auxiliary array z_aux ',       &
                 'failed; the program stops here (subroutine GeneratePart). '
               stop
            endif
         endif
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
         do i_vertex=Partz(Nz)%ID_first_vertex_sel,Partz(Nz)%ID_last_vertex_sel
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
            call z_min_max_DEM_DTM_9p_stencil(z_min_flag,Nz,i_vertex,          &
               z_aux(i_vertex))
! Generating daughter points of the Cartesian topography, according to dx
! (maybe the function "ceiling" would be better than "int" here)
               n_levels = int((h_reservoir(Nz) - z_aux(i_vertex)) / Domain%dx)
               do i=1,aux_factor
                  do j=1,aux_factor
                     do k=1,n_levels
                        call particle_position_extrusion(i_vertex,aux_factor,i,&
                           j,k,Partz(Nz)%Car_top_zone,h_reservoir(Nz),test_z,  &
                           pg_aux(NumParticles)%coord(1:3))
                        call particles_in_out_dams(Nz,                         &
                           pg_aux(NumParticles)%coord(1:3),test_dam)
! Update fluid particle counter. Check if the fluid particle is visible both to 
! the bottom reservoir and the optional dam.
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
!$omp shared(Domain,Nz,Mate,pg,Partz,max_rnd_eps)                              &
!$omp private(npi)
            do npi=Partz(Nz)%limit(1),Partz(Nz)%limit(2)
! To set particle parameters
               call SetParticleParameters(npi,Nz,Mate)
! To impose a white noise to particle positions
               if (Domain%RandomPos=='r') then
                  call pos_plus_white_noise(max_rnd_eps,pg(npi)%coord(1:3))
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
                  Xminreset(i) = MinOfMin(i) + NumPartPrima * Domain%dx
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
enddo
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
