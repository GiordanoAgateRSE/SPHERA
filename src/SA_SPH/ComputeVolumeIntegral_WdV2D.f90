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
! Program unit: ComputeVolumeIntegral_WdV2D                                    
! Description: Computing the integral of WdV extented to the volume delimited by the kernel support (radius=2h) of the particle i,
!              whose local coordinates are xpi=LocXY(1,icbs) and ypi=LocXY(2,icbs), and the adjacent boundary side icbs.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ComputeVolumeIntegral_WdV2D(icbs,Ncbslocal,Cloboside,LocXY,         &
   BoundarySide,xpmin,xpmax,interlen,VIntWdV)
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),parameter :: Nalfadiv = 10
double precision,parameter :: eps=0.05d0
integer(4),intent(IN) :: icbs,Ncbslocal
double precision,intent(IN) :: xpmin,xpmax,interlen
integer(4),intent(IN),dimension(1:MAXCLOSEBOUNDSIDES) :: Cloboside
double precision,intent(IN),dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY
type (TyBoundarySide),intent(IN),dimension(1:NumBSides) :: BoundarySide
double precision,intent(INOUT) :: VIntWdV
integer(4) :: jcbs,ndiv,ipt
double precision :: xpi,ypi,yplimite,ypj,angle,dalfarif,Intalfa,tanalfa,ris
double precision :: dalfa,alfaA,alfaB,alfa,csiPA,etaPA,csiPB,etaPB
character(len=lencard) :: nomsub = "ComputeVolumeIntegral_WdV2D"
double precision,external :: WIntegr
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
VIntWdV = zero
if (interlen<=zero) return
yplimite = eps * Domain%h
xpi = LocXY(1,icbs)
ypi = LocXY(2,icbs)
!------------------------
! Statements
!------------------------
if (ypi>=yplimite) then
   csiPA = xpmin - xpi
   etaPA = ypi
   csiPB = xpmax - xpi
   etaPB = ypi
   alfaA = Atan2(csiPA, etaPA)
   alfaB = Atan2(csiPB, etaPB)
   Intalfa = alfaB - alfaA
   ndiv = Int(intalfa / dalfarif + half)
   if (ndiv<2) ndiv = 2
   dalfa = intalfa / ndiv
   alfa = alfaA - half * dalfa
   do ipt=1,ndiv
      alfa = alfa + dalfa
      tanalfa = Tan(alfa)
      ris = ypi * Dsqrt(one + tanalfa * tanalfa)
      VIntWdV = VIntWdV + WIntegr(ris, Domain%h) * dalfa
   enddo
   elseif (ypi<yplimite.and.Ncbslocal==1) then        
! There is only one close side: "icbs"
      if (xpmin>=xpi.and.xpi<=xpmax) then      
         VIntWdV = half
         else    
            VIntWdV = zero
      endif
      elseif (ypi<yplimite.and.Ncbslocal==2) then
! There can be two sides close to the particle 
! Index of the second side, which is close to the particle 
         jcbs = Ncbslocal + 1 - icbs        
! Distance between the particle and the second side      
         ypj = LocXY(2, jcbs)                    
         if (ypj<=yplimite) then               
! The particle is very close to the vertex, which is in common between the sides 
            angle = zero
            if (BoundarySide(Cloboside(1))%previous_side==Cloboside(2)) then
               angle = BoundarySide(Cloboside(1))%angle
               elseif (BoundarySide(Cloboside(2))%previous_side==Cloboside(1)) &
                  then
                  angle = BoundarySide(Cloboside(2))%angle
                  else
                     write(nout,'(a,2i10)') 'ERROR!! Sides not consecutive',   &
                        Cloboside(1),Cloboside(2)
                     call diagnostic(arg1=8,arg2=4,arg3=nomsub)
            endif
            VIntWdV = half * (angle + PIGRECO) / (two * PIGRECO)
            else    
! The particle is very close only to the side "icbs"
               if (xpmin>=xpi.and.xpi<=xpmax) then      
                  VIntWdV = half
                  else    
                  VIntWdV = zero
               endif         
         endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeVolumeIntegral_WdV2D

