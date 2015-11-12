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
! Program unit: LocalNormalCoordinates  
! Description: Given the local coordinates PX(1 to 2) of a point P laying on the plane of the boundary face nf, the procedure
!              assigns to csi(1 to 3) the normal coordinates of the point Q corresponding to P in the inverse linear tranformation.  
!----------------------------------------------------------------------------------------------------------------------------------

subroutine LocalNormalCoordinates ( PX, csi, nf )
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
integer(4),intent(IN) :: nf
double precision,intent(IN), dimension(1:SPACEDIM) :: PX
double precision,intent(INOUT),dimension(1:SPACEDIM) :: csi
integer(4) :: i,j,k,nodes,fkod
double precision :: AA,BB,CC,DueArea,UsuDueArea,xj,yj,xk,yk
integer(4), dimension(3)   :: iseg = (/ 2,3,1 /)
integer(4), dimension(2,3) :: mainod 
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! Modification for compatibility xlf90 (IBM compiler): start
mainod(1,1) = 1
mainod(1,2) = 2
mainod(1,3) = 3
mainod(2,1) = 1
mainod(2,2) = 3
mainod(2,3) = 4
! Modification for compatibility xlf90 (IBM compiler): end
nodes = 4
if (BoundaryFace(nf)%Node(4)%name<=0) nodes = 3
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
end do
csi(3) = one - (csi(1) + csi(2))
!------------------------
! Deallocations
!------------------------
return
end subroutine LocalNormalCoordinates

