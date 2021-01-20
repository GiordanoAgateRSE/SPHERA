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
! Program unit: IsParticleInternal3D
! Description: To check whether a particle is internal to a 3D volume 
!              (domain/zone) or not. It checks if point "Px" is internal to  
!              the perimeter "mib". It returns ".true." (positive check) or 
!              ".false.". The perimeter can be both convex or concave. The 
!              input point is internal to the zone if the number of both "faces 
!              intercepted by the vertical (passing for the input point) above 
!              the input point" and "faces intercepted by the vertical (passing 
!              for the input point) below the input point" are odd.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
logical function IsParticleInternal3D(mib,PX,IsopraS)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),parameter :: intxy = 3
double precision,parameter :: eps = 0.001d0
integer(4),intent(in) :: mib
double precision,intent(in) :: PX(SPACEDIM)
integer(4),intent(in) :: IsopraS
integer(4) :: kf,nf,i,j,sd,nnodes,norig,Nints,IntSotto,IntSopra
integer(4) :: test
double precision :: tpar
double precision :: P1(SPACEDIM),Pint(SPACEDIM),LPint(SPACEDIM)
double precision,dimension(Tratto(mib)%numvertices) :: XYInts
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
   subroutine point_inout_quadrilateral(point,point_pol_1,point_pol_2,         &
                                        point_pol_3,point_pol_4,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_quadrilateral
   subroutine point_inout_pentagon(point,point_pol_1,point_pol_2,              &
                                   point_pol_3,point_pol_4,point_pol_5,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_pentagon
   subroutine point_inout_hexagon(point,point_pol_1,point_pol_2,               &
                                  point_pol_3,point_pol_4,point_pol_5,         &
                                  point_pol_6,test)
      implicit none
      double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2)
      double precision,intent(in) :: point_pol_3(2),point_pol_4(2)
      double precision,intent(in) :: point_pol_5(2),point_pol_6(2)
      integer(4),intent(inout) :: test
   end subroutine point_inout_hexagon
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
IsParticleInternal3D = .false.
Nints = 0
IntSotto = 0
IntSopra = 0
IntSopra = IsopraS
nnodes = 0
!------------------------
! Statements
!------------------------
! Loop on the zone faces: start
do kf=Tratto(mib)%iniface,(Tratto(mib)%iniface+Tratto(mib)%numvertices-1)
! Face ID
   nf = BFaceList(kf)
! Number of the vertices of the zone face   
   nnodes = 6
   if (BoundaryFace(nf)%Node(6)%name<=0) nnodes = 5
   if (BoundaryFace(nf)%Node(5)%name<=0) nnodes = 4
   if (BoundaryFace(nf)%Node(4)%name<=0) nnodes = 3
! The last vertex is arbitrarily assumed as the origin of the local reference 
! system 
   norig = nnodes 
   do sd=1,SPACEDIM
      P1(sd) = Vertice(sd,BoundaryFace(nf)%Node(norig)%name)
   enddo
! Relative distance between the input point (to test) and the zone face 
! (temporary value of tpar)
   tpar = zero
   do sd=1,SPACEDIM
      tpar = tpar + BoundaryFace(nf)%T(sd,3) * (P1(sd) - PX(sd))
   enddo
! To consider only non-vertical faces
   if (abs(BoundaryFace(nf)%T(3,3))>eps) then
! Relative distance between the input point (to test) and the zone face 
! along the vertical
      tpar = tpar / BoundaryFace(nf)%T(3,3)
! Pint: global coordinates of the intersection point between the zone face and 
!       the vertical line passing for the input point
      do sd=1,SPACEDIM 
         Pint(sd) = PX(sd)
      enddo
      Pint(3) = Pint(3) + tpar
! LPint: local coordinates of the intersection point between the zone face and 
!        the vertical line passing for the input point   
      LPint = zero
      do sd=1,PLANEDIM
         LPint(sd) = zero
         do j=1,SPACEDIM
            LPint(sd) = LPint(sd) + BoundaryFace(nf)%T(j,sd) * (Pint(j) - P1(j))
         enddo
      enddo
! If the face intercepts the vertical line passing for Px, it saves the 
! z-coordinate of the intersection point.
      test = 0
      select case (nnodes)
         case(3)
            call point_inout_convex_non_degenerate_polygon(LPint,nnodes,       &
               BoundaryFace(nf)%Node(1)%LX(1:2),                               &
               BoundaryFace(nf)%Node(2)%LX(1:2),                               &
               BoundaryFace(nf)%Node(3)%LX(1:2),                               &
               BoundaryFace(nf)%Node(3)%LX(1:2),                               &
               BoundaryFace(nf)%Node(3)%LX(1:2),                               &
               BoundaryFace(nf)%Node(3)%LX(1:2),test)         
         case(4)
            call point_inout_quadrilateral(LPint,                              &
                                           BoundaryFace(nf)%Node(1)%LX(1:2),   &
                                           BoundaryFace(nf)%Node(2)%LX(1:2),   &
                                           BoundaryFace(nf)%Node(3)%LX(1:2),   &
                                           BoundaryFace(nf)%Node(4)%LX(1:2),   &
                                           test)
         case(5)
            call point_inout_pentagon(LPint,BoundaryFace(nf)%Node(1)%LX(1:2),  &
                                      BoundaryFace(nf)%Node(2)%LX(1:2),        &
                                      BoundaryFace(nf)%Node(3)%LX(1:2),        &
                                      BoundaryFace(nf)%Node(4)%LX(1:2),        &
                                      BoundaryFace(nf)%Node(5)%LX(1:2),test)            
         case(6)
            call point_inout_hexagon(LPint,BoundaryFace(nf)%Node(1)%LX(1:2),   &
                                     BoundaryFace(nf)%Node(2)%LX(1:2),         &
                                     BoundaryFace(nf)%Node(3)%LX(1:2),         &
                                     BoundaryFace(nf)%Node(4)%LX(1:2),         &
                                     BoundaryFace(nf)%Node(5)%LX(1:2),         &
                                     BoundaryFace(nf)%Node(6)%LX(1:2),test)
         case default
            write(uerr,*) "Run-time error at IsParticleInternal3D. The number",&
               " of face vertices (nnodes) must be 3, 4, 5 or 6. SPHERA stops. "   
            stop                                      
      endselect
      if (test==1) then    
         Nints = Nints + 1
         XYInts(Nints) = Pint(3)
      endif
   endif
! Loop on the zone faces: end
enddo
! The input point is internal to the zone if the number of both "faces  
! intercepted by the vertical (passing for the input point) above the input   
! point" and "faces intercepted by the vertical (passing for the input point)  
! below the input point" are odd.
if (Nints>0) then
   do i=1,Nints
      if (XYInts(i)<=PX(intxy)) then
         IntSotto = IntSotto + 1
         else
            IntSopra = IntSopra + 1
      endif
   enddo
   if ((mod(IntSotto,2)==1).and.(mod(IntSopra,2)==1)) then
      IsParticleInternal3D = .true.
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end function IsParticleInternal3D
#endif
