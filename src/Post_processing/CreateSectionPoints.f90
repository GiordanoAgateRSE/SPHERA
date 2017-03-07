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
! Program unit: CreateSectionPoints               
! Description:      
!-------------------------------------------------------------------------------
subroutine CreateSectionPoints(vp,values,opt,seccor)
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
integer(4) :: seccor
double precision, dimension(3) :: vp
double precision, dimension(3,2) :: values
character(1) :: opt
integer(4) :: i,j,k,n,npse
integer(4),dimension(3) :: Nmesh
double precision,dimension(3) :: Cc,CcStart
character(100), external :: lcase
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! Creating a mesh of side dx
Nmesh = 1
npse = Control_Sections(seccor)%icont(1)
CcStart(:) = vp(:) - Domain%dx
!------------------------
! Statements
!------------------------
do n=1,SPACEDIM
   if ((n==1).AND.(lcase(opt)=="x")) cycle
   if ((n==2).AND.(lcase(opt)=="y")) cycle
   if ((n==3).AND.(lcase(opt)=="z")) cycle
   Nmesh(n) = nint((values(n,2) - values(n,1)) / Domain%dx )
   CcStart(n) = values(n,1) - Domain%dx * half
enddo
Cc(1) = CcStart(1)
do i=1,Nmesh(1)
   Cc(1) = Cc(1) + Domain%dx
   Cc(2) = CcStart(2)
   do j=1,Nmesh(2)
      Cc(2) = Cc(2) + Domain%dx
      Cc(3) = CcStart(3)
      do k=1,Nmesh(3)
         Cc(3) = Cc(3) + Domain%dx
         npse = npse + 1
         Control_Points(npse)%coord(:) = Cc(:)
      enddo
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine  CreateSectionPoints

