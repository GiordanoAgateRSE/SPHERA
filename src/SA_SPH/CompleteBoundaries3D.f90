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
! Program unit: CompleteBoundaries3D                                  
! Description: (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine CompleteBoundaries3D
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
integer(4) :: Kf,Nt,Nf
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
Kf = 0
BFaceList(1:NumFacce) = 0
!------------------------
! Statements
!------------------------
do Nt=1,NumTratti
   Tratto(Nt)%inivertex = 0
   Tratto(Nt)%numvertices = 0
   do Nf=1,NumFacce
      if (BoundaryFace(Nf)%stretch==Nt) then
         Kf = Kf + 1
         if (Tratto(Nt)%iniface==0) then 
! To store the initial face index
            Tratto(Nt)%iniface = Kf
         endif
         BFaceList(Kf) = Nf
         Tratto(Nt)%numvertices = Tratto(Nt)%numvertices + 1
      endif
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine CompleteBoundaries3D

