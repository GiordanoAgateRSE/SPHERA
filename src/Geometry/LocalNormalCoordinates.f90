!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: LocalNormalCoordinates  
! Description: Given the local Cartesian coordinates PX(1:2) of a point P laying on the plane of the boundary face nf, this
!              procedure assigns to csi(1:3) the non-Cartesian coordinates in a different local reference systems, whose axis are 
!              aligned with the face sides. Scaling provides coordinate values to go from 0 to 1 for internal points. This 
!              procedure applies for triangular faces and simplier works for rectangular faces (not for quadrilateral faces).   
!----------------------------------------------------------------------------------------------------------------------------------
subroutine LocalNormalCoordinates(PX,csi,nf)
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
integer(4),intent(in) :: nf
double precision,intent(in) :: PX(SPACEDIM)
double precision,intent(inout) :: csi(SPACEDIM)
integer(4) :: i,j,k,nodes,fkod
double precision :: AA,BB,CC,DueArea,UsuDueArea,xj,yj,xk,yk
integer(4),dimension(3) :: iseg = (/ 2,3,1 /)
integer(4) :: mainod(2,3)
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
mainod(1,1) = 1
mainod(1,2) = 2
mainod(1,3) = 3
mainod(2,1) = 1
mainod(2,2) = 3
mainod(2,3) = 4
nodes = 4
if (BoundaryFace(nf)%Node(4)%name<0) nodes = 3
fkod = nodes - 2                    ! = 1 (triangle), =2 (quadrilateral)
DueArea = (3 - fkod) * BoundaryFace(nf)%Area
UsuDueArea = one / DueArea
!------------------------
! Statements
!------------------------
do i=1,2
   j = iseg(i)
   k = iseg(j)
   xj = BoundaryFace(nf)%Node(mainod(fkod,j))%LX(1)
   yj = BoundaryFace(nf)%Node(mainod(fkod,j))%LX(2)
   xk = BoundaryFace(nf)%Node(mainod(fkod,k))%LX(1)
   yk = BoundaryFace(nf)%Node(mainod(fkod,k))%LX(2)
   AA = xj * yk - xk * yj
   BB = yj - yk
   CC = xk - xj
   csi(i) = (AA + BB * PX(1) + CC * PX(2)) * UsuDueArea
enddo
csi(3) = one - (csi(1) + csi(2))
!------------------------
! Deallocations
!------------------------
return
end subroutine LocalNormalCoordinates

