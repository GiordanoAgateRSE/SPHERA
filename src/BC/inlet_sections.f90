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
! Program unit: inlet_sections
! Description: Elaboration of the quantities of the inlet sections 
!-------------------------------------------------------------------------------
subroutine inlet_sections
!------------------------
! Modules
!------------------------
use Static_allocation_module
use I_O_file_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nt,isi,sd,ip,i_source,i_rec
#ifdef SPACE_3D
integer(4) :: i,j,NumPartR,NumPartS,nodes
#elif defined SPACE_2D
integer(4) :: nA
#endif
double precision :: deltapart,fluid_depth,weir_length
#ifdef SPACE_3D
double precision :: etalocal,eps,deltaR,deltaS,LenR,LenS,distR,distS,csi
#elif defined SPACE_2D
double precision :: linedist,sidelen,eps
double precision,dimension(1:SPACEDIM) :: A,ss
#endif
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! searching the ID of the source side
#ifdef SPACE_3D
SourceFace = 0
#elif defined SPACE_2D
SourceSide = 0
#endif
SpCount = 0
i_source=0
!------------------------
! Statements
!------------------------
#ifdef SPACE_3D
do isi=1,NumFacce
   nt = BoundaryFace(isi)%stretch
   if (Tratto(nt)%tipo=="sour") then
      SourceFace = isi
#elif defined SPACE_2D
do isi=1,NumBSides
   if (BoundarySide(isi)%tipo=="sour") then
      SourceSide = isi
#endif
      i_source = i_source + 1
#ifdef SPACE_3D
      if (SourceFace>0) then
         nt = BoundaryFace(SourceFace)%stretch
         mat = Tratto(nt)%Medium
         izone = Tratto(nt)%zone
! Note: to insert a check in case of nodes=/4
         nodes = BoundaryFace(SourceFace)%nodes
#elif defined SPACE_2D
      nt = BoundarySide(SourceSide)%stretch
      irz = Tratto(nt)%zone
      mat = partz(irz)%Medium 
      nA = BoundarySide(SourceSide)%Vertex(1)
      do sd=1,SPACEDIM
         A(sd) = Vertice(sd,nA)
         ss(sd) = BoundarySide(SourceSide)%T(sd,1)
         nn(sd) = BoundarySide(SourceSide)%T(sd,3)
      enddo
#endif
         deltapart = Domain%dx
#ifdef SPACE_3D
! LenR and LenS are the length scales of the inlet section: they are computed
! as the distance between the first and the last inlet vertices and the third
! and the last inlet vertices, respectively. Particles are aligned with Plast-P1
! and Plast-P3, where P1 the first boundary vertex, ..., Plast being the last 
! boundary vertex. In case of a triangular inlet, we have particles aligned 
! with one direction: P3-P1. In case of a quadrilateral inlet, we have particles 
! distributed along two directions: P4-P1 and P4-P3.
         LenR = zero
         LenS = zero
         do sd=1,SPACEDIM
            LenR = LenR + (BoundaryFace(SourceFace)%Node(1)%GX(sd) -           &
                   BoundaryFace(SourceFace)%Node(nodes)%GX(sd)) ** 2
            LenS = LenS + (BoundaryFace(SourceFace)%Node(3)%GX(sd) -           &
                   BoundaryFace(SourceFace)%Node(nodes)%GX(sd)) ** 2
         enddo
         LenR = dsqrt(LenR)
         LenS = dsqrt(LenS)
#endif
         if (Tratto(nt)%time_flag.eqv..true.) then
! Linear time interpolation for the inlet flow rate: start
            do i_rec=1,Tratto(nt)%n_time_records
               if (Tratto(nt)%time_records(i_rec,1)>=simulation_time) then
                  if (Tratto(nt)%time_records(i_rec,1)==simulation_time) then
                     Tratto(nt)%FlowRate = Tratto(nt)%time_records(i_rec,2)
                     else
                        Tratto(nt)%FlowRate =                                  &
                           Tratto(nt)%time_records(i_rec-1,2)                  &
                           + (Tratto(nt)%time_records(i_rec,2) -               &
                           Tratto(nt)%time_records(i_rec-1,2)) /               &
                           (Tratto(nt)%time_records(i_rec,1) -                 &
                           Tratto(nt)%time_records(i_rec-1,1)) *               &
                           (simulation_time -                                  &
                           Tratto(nt)%time_records(i_rec-1,1))
                  endif
                  exit
               endif
            enddo
! Linear time interpolation for the inlet flow rate: end
! Linear time interpolation for the inlet fluid depth: start
            if (Tratto(nt)%weir_flag.eqv..false.) then
               do i_rec=1,Tratto(nt)%n_time_records
                  if (Tratto(nt)%time_records(i_rec,1)>=simulation_time) then
                     if (Tratto(nt)%time_records(i_rec,1)==simulation_time) then
                        fluid_depth = Tratto(nt)%time_records(i_rec,3)
                        else
                           fluid_depth = Tratto(nt)%time_records(i_rec-1,3)    &
                              + (Tratto(nt)%time_records(i_rec,3) -            &
                              Tratto(nt)%time_records(i_rec-1,3)) /            &
                              (Tratto(nt)%time_records(i_rec,1) -              &
                              Tratto(nt)%time_records(i_rec-1,1)) *            &
                              (simulation_time -                               &
                              Tratto(nt)%time_records(i_rec-1,1))
                     endif
                     exit
                  endif
               enddo
               else            
! LRRA approximation of Swamee (1988, JHE) formula for Bazin weirs with 
! upstream velocity
#ifdef SPACE_3D
                  weir_length = LenS
#elif defined SPACE_2D
               weir_length = 1.d0
#endif
                  fluid_depth = 3.d0 * ((Tratto(nt)%FlowRate / weir_length) ** &
                     (2.d0/3.d0)) / (9.806d0 ** (1.d0/3.d0))
            endif
#ifdef SPACE_3D
            BoundaryFace(SourceFace)%Node(1:2)%GX(3) =                         &
               BoundaryFace(SourceFace)%Node(3)%GX(3) + fluid_depth
#elif defined SPACE_2D
         BoundarySide(SourceSide)%length = fluid_depth
#endif
! Linear time interpolation for the inlet fluid depth: end
         endif      
#ifdef SPACE_3D   
         NumPartR = int(LenR / deltapart + 0.01d0)
         NumPartS = int(LenS / deltapart + 0.01d0)
         deltaR = LenR / NumPartR 
         deltaS = LenS / NumPartS 
         eps = -half
         zfila = eps * deltapart
         distR = -half * deltaR
         ip = 0
         do i=1,NumPartR
            distR = distR + deltaR
            csi = distR / LenR
            distS = -half * deltaS
            do j=1,NumPartS
               distS = distS + deltaS
               etalocal = distS / LenS
               ip = ip + 1
               do sd=1,SPACEDIM
                  P(sd) = BoundaryFace(SourceFace)%Node(4)%GX(sd) * (one -     &
                          csi) + BoundaryFace(SourceFace)%Node(1)%GX(sd) * csi
                  Q(sd) = BoundaryFace(SourceFace)%Node(3)%GX(sd) * (one -     &
                          csi) + BoundaryFace(SourceFace)%Node(2)%GX(sd) * csi
                  PartLine(i_source,ip,sd) = P(sd) * (one - etalocal) + Q(sd) *&
                                             etalocal
               enddo
            enddo
         enddo
         NumPartFace(i_source) = ip
#elif defined SPACE_2D
      sidelen = BoundarySide(SourceSide)%length
      NumPartperLine(i_source) = int(sidelen / deltapart + 0.01d0)
      eps = -half
      yfila = eps * deltapart
      linedist = -half * deltapart
      do ip=1,NumPartperLine(i_source)
         linedist = linedist + deltapart
         do sd=1,SPACEDIM
            PartLine(i_source,ip,sd) = A(sd) + linedist * ss(sd)
         enddo
      enddo
#endif
         ParticleVolume = Domain%PVolume
#ifdef SPACE_3D
         RowPeriod = ParticleVolume * NumPartFace(i_source) /                  &
                     Tratto(nt)%FlowRate
#elif defined SPACE_2D
      RowPeriod = ParticleVolume * NumPartperLine(i_source) /                  &
                  Tratto(nt)%FlowRate 
#endif
         RowVelocity(i_source) = Domain%dx / RowPeriod
         Tratto(nt)%NormVelocity = RowVelocity(i_source)
#ifdef SPACE_2D
      partz(irz)%vel(1) = RowVelocity(i_source) * nn(1)
      partz(irz)%vel(2) = RowVelocity(i_source) * nn(2)
      partz(irz)%vel(3) = RowVelocity(i_source) * nn(3)
#endif
         pinttimeratio = -1
#ifdef SPACE_3D
      endif
#endif
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine inlet_sections
