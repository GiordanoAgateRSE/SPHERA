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
! Program unit: ComputeSurfaceIntegral_WdS2D                                   
! Description:  Computing the surface integral of kernel W along the segments intercepted by the kernel support (radius=2h) of the
!               particle i, whose local coordinates are xpi=LocXY(1,icbs) and ypi=LocXY(2,icbs), on the adjacent boundary 
!               side icbs.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ComputeSurfaceIntegral_WdS2D(icbs,LocXY,xpmin,interlen,SIntWds)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4), parameter :: ndivrif = 2
integer(4),intent(IN) :: icbs
double precision,intent(IN) :: xpmin,interlen
double precision,intent(IN),dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
double precision,intent(INOUT) :: SIntWds
integer(4) :: ndiv,ipt 
double precision :: xpi,ypi,ypiq,dsrif,deltas,xpip,dxpip,ris,risq,Wds
double precision,external :: w
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
SIntWds = zero
dsrif = Domain%h / ndivrif
!------------------------
! Statements
!------------------------
if (interlen>zero) then
   xpi = LocXY(1, icbs)
   ypi = LocXY(2, icbs)
   if (ypi<zero) then      
      ypi = zero 
   endif    
   ypiq = ypi * ypi
   ndiv = Int(interlen / dsrif + half)
   if (ndiv<ndivrif) then
      ndiv = ndivrif
   endif
   deltas = interlen / ndiv
   xpip = xpmin - half * deltas
   do ipt=1,ndiv
      xpip = xpip + deltas
      dxpip = xpip - xpi
      risq = dxpip * dxpip + ypiq
      ris = Dsqrt(risq)
      Wds = w(ris, Domain%h, Domain%coefke) * deltas
      SIntWds = SIntWds + Wds
   enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeSurfaceIntegral_WdS2D

