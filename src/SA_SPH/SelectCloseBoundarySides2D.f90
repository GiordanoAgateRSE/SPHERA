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
! Program unit: SelectCloseBoundarySides2D                                          
! Description:  Selecting among the close boundary sides, those that really give contribution to the equations of particle 'npi'.
!               It returns:
!                  IntNcbs: number of close boundary sides, which give contribution (= 0, 1, 2)
!                  Intboside(1:IntNcbs): list of close boundary sides, which give contribution
!                  IntLocXY(1:PLANEDIM,1:Ncbs): local coordinates of particle np with respect each boundary side, which gives
!                                               contribution                                
!               (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine SelectCloseBoundarySides2D(npi,Ncbs,Cloboside,LocXY,IntNcbs,        &
   Intboside,IntLocXY)
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
double precision, parameter :: eps = 0.501d0
integer(4),intent(in) :: npi,Ncbs
integer(4),dimension(1:MAXCLOSEBOUNDSIDES),intent(in) :: Cloboside
double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES),intent(in)  :: LocXY
integer(4),intent(out) :: IntNcbs
integer(4),dimension(1:MAXCLOSEBOUNDSIDES),intent(out) :: Intboside
double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES),intent(out) :: IntLocXY
integer(4) :: icbs,isi,nt
double precision :: sidelen,yxpmin
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
IntNcbs = 0
Intboside(:) = 0
IntLocXY(:,:) = zero
!------------------------
! Statements
!------------------------
! Loop over the close boundaries previously found
do icbs=1,Ncbs
   isi = Cloboside(icbs)
   sidelen = BoundarySide(isi)%length
   nt = BoundarySide(isi)%stretch
! The boundary is of type "source": the side is considered only if the particle 
! is in front of the side itself and inside the domain or very closest the 
! side (<0.5*Domain%dd)
   if (Tratto(nt)%tipo=="sour") then
      yxpmin = -eps * Domain%dd
      if ((LocXY(2,icbs)>yxpmin).and.(LocXY(1,icbs)>zero).and.                 &
         (LocXY(1,icbs)<sidelen)) then
         IntNcbs = IntNcbs + 1
         Intboside(IntNcbs) = Cloboside(icbs)
         IntLocXY(:,IntNcbs) = LocXY(:,icbs)
      endif
! The boundary is of type "outlet velocity"
      elseif ((Tratto(nt)%tipo=="velo").or.(Tratto(nt)%tipo=="flow")) then
         yxpmin = zero
         if ((LocXY(2,icbs)>yxpmin).and.(LocXY(1,icbs)>zero).and.              &
            (LocXY(1,icbs)<sidelen)) then
            IntNcbs = IntNcbs + 1
            Intboside(IntNcbs) = Cloboside(icbs)
            IntLocXY(:,IntNcbs) = LocXY(:,icbs)
! On-going particle close to an outlet section (see erosion criterion)
            pg(npi)%CloseBcOut = 1
         endif
! The boundary is of type "open outlet"
         elseif (Tratto(nt)%tipo=="open") then
            yxpmin = zero
            if ((LocXY(2,icbs)>yxpmin).and.(LocXY(1,icbs)>zero).and.           &
               (LocXY(1,icbs)<sidelen)) then
               IntNcbs = IntNcbs + 1
               Intboside(IntNcbs) = Cloboside(icbs)
               IntLocXY(:,IntNcbs) = LocXY(:,icbs)
! On-going particle close to an outlet section (see erosion criterion)
               pg(npi)%CloseBcOut = 1
            endif
! The boundary type is neither "source", nor "velo", "flow", "open"
            else
               yxpmin = zero          
! On-going particle is inside the domain and in front of a boundary of 
! type "open"
               if (LocXY(2,icbs)>yxpmin) then
! The closest boundary is accounted  
                  IntNcbs = IntNcbs + 1
                  Intboside(IntNcbs) = Cloboside(icbs)
                  IntLocXY(:,IntNcbs) = LocXY(:,icbs)
               endif
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine SelectCloseBoundarySides2D

