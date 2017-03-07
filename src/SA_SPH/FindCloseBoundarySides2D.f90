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
! Program unit: FindCloseBoundarySides2D                                       
! Description: To find the "close" boundary sides, i.e. those sited at a 
!              distance from particle npi<=2h. It returns:
!                 Ncbs: number of close boundary sides (= 0, 1, 2)
!                 Cloboside(1:Ncbs): list of close boundary sides
!                 LocXY(1:PLANEDIM,1:Ncbs): local coordinates of particle npi 
!                                           with respect each boundary side 
!                                           (vertex V1) 
!              (Di Monaco et al., 2011, EACFM)
!-------------------------------------------------------------------------------
subroutine FindCloseBoundarySides2D(npi,Ncbs,Cloboside,LocXY)
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
integer(4),intent(out) :: Ncbs
integer(4),dimension(1:MAXCLOSEBOUNDSIDES),intent(out) :: Cloboside
double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES),intent(out) :: LocXY
integer(4) :: icbs,icb,isous,isi,v1,v2,sd,iside1,iside2 
double precision :: xp,yp,sidel,xpmin,xpmax,ypmin,ypmax,xpq
double precision :: Lmxpq,sidelen,ypmn,ypmx
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: Plocal,P1,P1P,sss,nnn
character(len=lencard) :: nomsub = "FindCloseBoundarySides2D"
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
acix(1) = 1  
acix(2) = 3
Cloboside = 0
LocXY = zero
Ncbs  = 0
Plocal(:) = pg(npi)%coord(acix(:))
ypmin = - doubleh
ypmax = doubleh
pg(npi)%CloseBcOut = 0
!------------------------
! Statements
!------------------------
! Loop over all the boundary sides of the domain
side_loop: do isi=1,NumBSides
! The boundary types "perimeter" and "pool" are skipped
   if (BoundarySide(isi)%tipo/="peri".AND.BoundarySide(isi)%tipo/="pool") then
! To load the side coordinates P1 and the local directional cosines "sss" and 
! "nnn"
! To estimate the distances between the current particle position "Plocal"  
! and the reference vertex "V1"
      v1 = BoundarySide(isi)%Vertex(1)
      v2 = BoundarySide(isi)%Vertex(2)
      do sd=1,PLANEDIM
         P1(sd) = Vertice(acix(sd),v1)
         P1P(sd) = Plocal(sd) - P1(sd)
         sss(sd) = BoundarySide(isi)%T(acix(sd), acix(1))
         nnn(sd) = BoundarySide(isi)%T(acix(sd), acix(2))
      enddo
! Particle coordinates "xp" and "yp" in the side local reference 
! system, where x-axis is aligned with the side segment and y-axis is normal 
! to the vector of vertices "V1" and the origin.
      sidel = BoundarySide(isi)%length
      xp = P1P(1) * sss(1) + P1P(2) * sss(2)
      yp = P1P(1) * nnn(1) + P1P(2) * nnn(2)
! To set the interaction area in the neighbourhood of the side segment having a 
! distance equal to +/- 2h in the y local direction and a distance of 2h from 
! the vertices V1 and V2 in the x local direction
      xpmin = - doubleh
      xpmax = sidel + doubleh
! To check if the particle has local coordinates falling inside the 
! interaction area
      if (((xpmin<xp).AND.(xp<xpmax)).AND.((ypmin<yp).AND.(yp<ypmax))) then
         if (xp<zero) then
            xpq = xp * xp
            ypmx = Dsqrt(doublesquareh - xpq)
            ypmn = -ypmx 
! The particle falls on the segment 
            elseif (xp<=sidel) then
               ypmx = ypmax    
               ypmn = ypmin 
               elseif (xp<xpmax) then
                  Lmxpq = (sidel - xp) * (sidel - xp)
                  ypmx = Dsqrt(doublesquareh - Lmxpq)
                  ypmn = -ypmx 
         endif
! The boundary must be considered for the current particle as a close boundary
         if ((ypmn<yp).AND.(yp<ypmx)) then
! The number of close boundaries is increased
            Ncbs = Ncbs + 1
! To check the maximum number allowed for closest boundaries
            if (Ncbs>MAXCLOSEBOUNDSIDES) then
! The particle "npi" has more than two boundary sides: to reduce dx and restart
               write(nout,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,       &
                  ' has more than two boundary sides: reduce dd and restart.'
               write(nout,'(1x,a,3f15.10)') '        Coordinate: ',            &
                  pg(npi)%coord(:)
               write(nscr,'(1x,a,i12,a)')   ' ERROR! The particle ',npi,       &
                  ' has more than two boundary sides: reduce dd and restart.'
               write(nscr,'(1x,a,3f15.10)') '        Coordinate: ',            &
                  pg(npi)%coord(:)
               call diagnostic(arg1=8,arg2=8,arg3=nomsub)
            endif
! To save the boundary ID and the local coordinates on the boundary segment
            Cloboside(Ncbs) = isi
            LocXY(1,Ncbs) = xp
            LocXY(2,Ncbs) = yp
         endif
      endif
   endif
enddo side_loop
! Searching a nearest side of type "source" 
isous = 0
do icbs=1,Ncbs
   isi = Cloboside(icbs)
! To increase the number of particles close to boundaries and to compute the 
! maximum heigth 
!$omp critical (numpa)
   BoundarySide(isi)%CloseParticles = BoundarySide(isi)%CloseParticles + 1
   if (BoundarySide(isi)%CloseParticles_maxQuota<pg(npi)%coord(3))             &
      BoundarySide(isi)%CloseParticles_maxQuota = pg(npi)%coord(3)
!$omp end critical (numpa)
   if (BoundarySide(isi)%tipo=="sour") then
      XP = LocXY(1,icbs)
      sidel = BoundarySide(isi)%length
      if ((XP>zero).AND.(XP<sidel)) then
! To mark the inlet section to erase the other possible close sides 
         isous = icbs    
         exit
         else                
! To cancel the inlet section, which has no influence here 
            if (icbs<Ncbs) then
               do icb=(icbs+1),Ncbs
                  Cloboside(icb-1) = Cloboside(icb)
                  LocXY(1,icb-1) = LocXY(1,icb)
                  LocXY(2,icb-1) = LocXY(2,icb)
               enddo
               Ncbs = Ncbs -1
               exit
               elseif (icbs==Ncbs) then
                  Ncbs = Ncbs -1
                  exit
            endif
      endif
   endif
enddo
! An inlet section ("source") has been found: the other nearest sides are 
! deleted
if (isous>0) then   
   Ncbs = 1
   Cloboside(Ncbs) = Cloboside(isous)
   LocXY(1,Ncbs) = LocXY(1,isous)
   LocXY(2,Ncbs) = LocXY(2,isous)
endif
! To check if more than two sides are close to the particle
if (Ncbs>LIMCLOSEBOUNDSIDES) then
   call diagnostic(arg1=8,arg2=9,arg3=nomsub)
endif
! There are two close boundaries: to change the origin of the local abscissa in 
! one of the two adjacent boundary sides, so that the abscissa origin coincides 
! with the common vertex of both in each side
if (Ncbs==2) then
   iside1 = Cloboside(1)
   iside2 = Cloboside(2)
   if (BoundarySide(iside1)%previous_side==iside2) then
      sidelen = BoundarySide(iside2)%length
      LocXY(1,2) = sidelen - LocXY(1,2)
      elseif (BoundarySide(iside2)%previous_side==iside1) then
         sidelen = BoundarySide(iside1)%length
         LocXY(1,1) = sidelen - LocXY(1,1)
         else
! Sides are not adjacent: this case should never occur if the minimum length
! of boundary sides is bigger than 4h 
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine FindCloseBoundarySides2D

