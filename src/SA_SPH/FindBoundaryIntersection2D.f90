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
! Program unit: FindBoundaryIntersection2D                                       
! Description: To find the intersection segment between the kernel support of particle i, whose local coordinates are
!              xpi=LocXY(1,icbs) and ypi=LocXY(2,icbs), and the straight boundary side iside=Cloboside(icbs), which lies on the 
!              local x-axis and extends from x=0 to bsidelen = BoundarySide(iside)%Length. It returns:
!                 xpmin: minimum abscissa of intersected segment
!                 xpmax: maximum abscissa of intersected segment
!                 interlen: length of the intersected segment   
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine FindBoundaryIntersection2D(icbs,Cloboside,LocXY,BoundarySide,xpmin, &
   xpmax,interlen)
   !------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: icbs
integer(4),dimension(1:MAXCLOSEBOUNDSIDES) :: Cloboside
double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
type (TyBoundarySide),dimension(1:NumBSides) :: BoundarySide
integer(4) :: iside
double precision :: xpmin,xpmax,interlen
double precision :: xpi,ypi,ypiq,bsidelen,halfchord,minlen,eps,XA,XB
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
eps = 0.01d0
minlen = eps*Domain%h
interlen = zero
xpi = LocXY(1,icbs)
ypi = LocXY(2,icbs)
if (ypi<zero) then    
! If the particle is out of the fluid domain, it is considered as if 
! it were on the boundary
   ypi = zero            
endif    
!------------------------
! Statements
!------------------------
ypiq = ypi * ypi
halfchord = Dsqrt(doublesquareh - ypiq)
XA = xpi - halfchord
XB = xpi + halfchord
iside = Cloboside(icbs)
bsidelen = BoundarySide(iside)%Length
xpmin = zero
if (XA>xpmin) then
   xpmin = XA
endif
xpmax = bsidelen
if (XB<xpmax) then
   xpmax = XB
endif
interlen = xpmax - xpmin
if (interlen<=minlen) then
! Intersection is too small
   interlen = zero               
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine FindBoundaryIntersection2D

