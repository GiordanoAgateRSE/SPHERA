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
! Program unit: BoundaryVolumeIntegrals2D                                 
! Description: To compute the boundary volume integrals IntWdV. 
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine BoundaryVolumeIntegrals2D(icbs,LocXY,xpmin,xpmax,interlen,IntWdV,   &
   IntDpWdV)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),parameter :: Nalfadiv = 10
double precision,parameter :: eps = 0.05d0
integer(4),intent(IN) :: icbs
double precision,intent(IN) :: xpmin,xpmax,interlen
double precision,intent(IN),dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
double precision,intent(INOUT) :: IntWdV
double precision,intent(INOUT),dimension(1:2) :: IntDpWdV
integer(4) :: k,ndiv
double precision :: xpi,ypi,yplimite,dalfarif,Intalfa,dalfa,alfaA,alfaB,alfa_k
double precision :: csiPA,csiPB,sinalfa,cosalfa,rb,rob,mult
double precision,external :: WIntegr,J2Wro2
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
dalfarif = PIGRECO / Nalfadiv
IntWdV = zero
IntDpWdV = zero
!------------------------
! Statements
!------------------------
if (interlen<=zero) return
yplimite = eps * Domain%h
xpi = LocXY(1, icbs)
ypi = LocXY(2, icbs)
if (ypi<yplimite) ypi = yplimite
csiPA = xpmin - xpi
csiPB = xpmax - xpi
alfaA = Atan2(csiPA, ypi)
alfaB = Atan2(csiPB, ypi)
Intalfa = alfaB - alfaA
ndiv = Int(intalfa / dalfarif + half)
if (ndiv<1) ndiv = 1
dalfa = intalfa / ndiv
alfa_k = alfaA - half * dalfa
do k=1,ndiv
   alfa_k = alfa_k + dalfa
   sinalfa = sin(alfa_k)
   cosalfa = cos(alfa_k)
   rb = ypi / cosalfa
   rob = rb / Domain%h
   mult = Domain%h * J2Wro2(rob) * dalfa
   IntDpWdV(1) = IntDpWdV(1) + sinalfa * mult
   IntDpWdV(2) = IntDpWdV(2) - cosalfa * mult
   IntWdV = IntWdV + WIntegr(rb, Domain%h) * dalfa
enddo
! To check if the particle is very close the the internal side  
if (ypi==yplimite) then
   IntWdV = half
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine BoundaryVolumeIntegrals2D

