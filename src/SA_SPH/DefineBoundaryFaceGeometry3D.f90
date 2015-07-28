!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: DefineBoundaryFaceGeometry3D                                    
! Description: To define boundary faces from 3D geometry.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine DefineBoundaryFaceGeometry3D
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
integer(4) :: Kf,Nf,Nt
double precision,dimension(SPACEDIM) :: Fi
double precision,dimension(SPACEDIM,SPACEDIM) :: TT,RGP,RMF
Data Fi /1.0d0,1.0d0,1.0d0/
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
do Kf=1,NumFacce,1
   Nf = BFaceList(Kf)
   if (Nf==0) cycle
   Nt = BoundaryFace(Nf)%stretch
   call DefineLocalSystemVersors (Nf)
   TT(1:SPACEDIM,1:SPACEDIM) = BoundaryFace(nf)%T(1:SPACEDIM,1:SPACEDIM)
   RGP = zero
   RMF = zero
   call BoundaryMassForceMatrix3D (TT, RMF, Fi)
   BoundaryFace(nf)%RFi(1:SPACEDIM,1:SPACEDIM) = RMF(1:SPACEDIM,1:SPACEDIM)
   if (Tratto(nt)%tipo=="tapi") then
      BoundaryFace(nf)%velocity(1:SPACEDIM) = Tratto(nt)%velocity(1:SPACEDIM)
      else
         BoundaryFace(nf)%velocity(1:SPACEDIM) = zero
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine DefineBoundaryFaceGeometry3D

