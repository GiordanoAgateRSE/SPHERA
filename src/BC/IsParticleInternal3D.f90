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
! Program unit: IsParticleInternal3D
! Description: To check whether a particle is internal to the 3D domain or not. It checks if point Px() is internal to the 
!              perimeter mib. It returns 'true' (positive check) or 'false'. The perimeter can be both convex or concave. 
!----------------------------------------------------------------------------------------------------------------------------------

Logical Function IsParticleInternal3D (mib,PX,IsopraS)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),parameter :: intxy = 3
double precision,parameter :: eps = 0.001d0
integer(4), intent(IN) :: mib
double precision,intent(IN),dimension(SPACEDIM) :: PX
integer(4), intent(IN) :: IsopraS
integer(4) :: kf,nf,i,j,sd,nnodes,norig,Nints,IntSotto,IntSopra,fkod
double precision :: tpar
double precision,dimension(SPACEDIM) :: P1, Pint, LPint
double precision,dimension(3) :: csi
double precision,dimension(Tratto(mib)%numvertices) :: XYInts
logical, external :: IsPointInternal
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
Nints = 0
IntSotto = 0
IntSopra = 0
IntSopra = IsopraS
IsParticleInternal3D = .FALSE.
!------------------------
! Statements
!------------------------
do kf=Tratto(mib)%iniface,(Tratto(mib)%iniface+Tratto(mib)%numvertices-1)
   nf = BFaceList(kf)
   nnodes = 4
   if (BoundaryFace(nf)%Node(4)%name<=0) nnodes = 3
! Node, which is the origin of the local system 
   norig = nnodes 
   do sd=1,SPACEDIM
      P1(sd) = Vertice(sd,BoundaryFace(nf)%Node(norig)%name)
   end do
   tpar = zero
   do sd=1,SPACEDIM
      tpar = tpar + BoundaryFace(nf)%T(sd,3) * (P1(sd) - PX(sd))
   end do
   if (Abs(BoundaryFace(nf)%T(3,3))>eps) then
      tpar = tpar / BoundaryFace(nf)%T(3,3)
! Pint: global coordinates of the intersection point
      do sd=1,SPACEDIM 
         Pint(sd) = PX(sd)
      end do
      Pint(3) = Pint(3) + tpar
      LPint = zero
! LPint: local coordinates of the intersection point
      do sd=1,PLANEDIM
         LPint(sd) = zero
         do j=1,SPACEDIM
            LPint(sd) = LPint(sd) + BoundaryFace(nf)%T(j,sd) * (Pint(j) -      &
                        P1(j))
         end do
      end do
      call LocalNormalCoordinates(LPint,csi,nf)
      fkod = nnodes - 2
! The intersection point is internal to the face; i.e. the face intersepts 
! the vertical line passing for Px. It saves the coordinate z of intersection.
      if (IsPointInternal(fkod,csi)) then    
         Nints = Nints + 1
         XYInts(Nints) = Pint(3)
      end if
   end if
end do
if (Nints>0) then
   do i=1,Nints
      if (XYInts(i)<=PX(intxy)) then
         IntSotto = IntSotto + 1
         else
            IntSopra = IntSopra + 1
      end if
   end do
   if ((Mod(IntSotto,2)==1).AND.(Mod(IntSopra,2)==1)) then
      IsParticleInternal3D = .TRUE.
   end if
end if
!------------------------
! Deallocations
!------------------------
return
end Function IsParticleInternal3D

