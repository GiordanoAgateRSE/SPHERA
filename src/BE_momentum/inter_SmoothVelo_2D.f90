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
! Program unit: inter_SmoothVelo_2D 
! Description: To calculate a corrective term for velocity.    
!-------------------------------------------------------------------------------
subroutine inter_SmoothVelo_2D 
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
integer(4) :: npi,i,ii,npj,contj,npartint,icbs,Ncbs,IntNcbs,ibdt,ibdp,iside
integer(4) :: sidestr
double precision :: rhoi,rhoj,amassj,pesoj,moddervel,unity,IntWdV   
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:PLANEDIM) :: sss,nnn,DVLoc,DVGlo,BCLoc,BCGlo
double precision,dimension(1:PLANEDIM) :: IntLocXY
double precision,dimension(3) :: dervel     
double precision,dimension(:),allocatable :: unity_vec
double precision,dimension(:,:),allocatable :: dervel_mat
character(4) :: strtype
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
if (n_bodies>0) then  
   allocate(dervel_mat(nag,3))
   allocate(unity_vec(nag))
   dervel_mat = 0.
   unity_vec = 0.
endif
acix(1) = 1
acix(2) = 3
unity = zero
!------------------------
! Statements
!------------------------
! Body particle contributions to pressure smoothing
if (n_bodies>0) then
   call start_and_stop(3,7)
   call start_and_stop(2,19)
   call body_to_smoothing_vel(dervel_mat,unity_vec)
   call start_and_stop(3,19)
   call start_and_stop(2,7)
endif
!$omp parallel do default(none)                                                &
!$omp private(ii,npi,contj,npartint,npj,rhoi,rhoj,amassj,dervel,moddervel)     &
!$omp private(pesoj,Ncbs,IntNcbs,ibdt,ibdp,icbs,IntLocXY,iside,sidestr,strtype)&
!$omp private(i,sss,nnn,DVGlo,DVLoc,IntWdV,BCLoc,BCGlo)                        &
!$omp shared(nag,pg,Domain,Med,Tratto,acix,nPartIntorno,NMAXPARTJ,PartIntorno) &
!$omp shared(PartKernel,BoundaryDataPointer,BoundaryDataTab,BoundarySide)      &
!$omp shared(indarrayFlu,Array_Flu,esplosione,kernel_fw,unity,dervel_mat)      &
!$omp shared(unity_vec,n_bodies,Granular_flows_options)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   pg(npi)%var = zero
! The mixture particles, which are temporarily affected by the frictional 
! viscosity threshold are fixed.
   if (pg(npi)%mu==Med(pg(npi)%imed)%mumx) cycle
   pg(npi)%Envar = zero
   do contj=1,nPartIntorno(npi)
      npartint = (npi-1)* NMAXPARTJ + contj
      npj = PartIntorno(npartint)
      rhoi = pg(npi)%dens
      rhoj = pg(npj)%dens
      amassj = pg(npj)%mass
      dervel(:) = pg(npj)%vel(:) - pg(npi)%vel(:)
      if (pg(npj)%vel_type/="std") then  
         rhoj = rhoi
         amassj = pg(npi)%mass
         moddervel = - two * (pg(npi)%vel(1) * pg(npj)%zer(1) + pg(npi)%vel(2) &
                     * pg(npj)%zer(2) + pg(npi)%vel(3) * pg(npj)%zer(3))
         dervel(:) = moddervel * pg(npj)%zer(:) 
      endif
      if (Med(pg(npj)%imed)%den0/=Med(pg(npi)%imed)%den0) cycle
      pesoj = amassj * PartKernel(4,npartint) / rhoj
      unity = unity + pesoj
      pg(npi)%var(:) = pg(npi)%var(:) + dervel(:) * pesoj   
      if (esplosione) pg(npi)%Envar = pg(npi)%Envar + (pg(npj)%IntEn -         &
                                      pg(npi)%IntEn) * pesoj
   enddo
   if (n_bodies>0) then
      pg(npi)%var(:) = pg(npi)%var(:) + dervel_mat(npi,:)
      unity = unity + unity_vec(npi)
   endif
   if (Domain%tipo=="bsph") then
      pg(npi)%var(:) = pg(npi)%var(:)
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
      call DBSPH_inlet_outlet(npi)
      else
         Ncbs = BoundaryDataPointer(1,npi)
         IntNcbs = BoundaryDataPointer(2,npi)
         ibdt = BoundaryDataPointer(3,npi)
         if (IntNcbs>0) then
            do icbs=1,IntNcbs
               ibdp = ibdt + icbs - 1
               IntLocXY(1:PLANEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
               iside = BoundaryDataTab(ibdp)%CloBoNum
               sidestr = BoundarySide(iside)%stretch
               strtype = Tratto(sidestr)%tipo
               if ((strtype=='sour').or.(strtype=='velo').or.(strtype=='flow'))&
                  then
                  pg(npi)%var(:) = zero   
                  exit  
               endif
               do i=1,PLANEDIM
                  sss(i) = BoundarySide(iside)%T(acix(i),1)
                  nnn(i) = BoundarySide(iside)%T(acix(i),3)
                  DVGlo(i) = two * (Tratto(sidestr)%velocity(acix(i)) -        &
                             pg(npi)%vel(acix(i)))
               enddo
               DVLoc(1) = sss(1) * DVGlo(1) + sss(2) * DVGlo(2)
               DVLoc(2) = nnn(1) * DVGlo(1) + nnn(2) * DVGlo(2)
               IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
               if ((strtype=='fixe').or.(strtype=='tapi')) then
                  BCLoc(1) = DVLoc(1) * IntWdV * Tratto(sidestr)%ShearCoeff
                  BCLoc(2) = DVLoc(2) * IntWdV
                  BCGlo(1) = sss(1) * BCLoc(1) + nnn(1) * BCLoc(2)
                  BCGlo(2) = sss(2) * BCLoc(1) + nnn(2) * BCLoc(2)
                  pg(npi)%var(1) = pg(npi)%var(1) + BCGlo(1)   
                  pg(npi)%var(3) = pg(npi)%var(3) + BCGlo(2)   
               endif
            enddo
         endif
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
if (n_bodies>0) then
   deallocate(dervel_mat)
   deallocate(unity_vec)
endif
return
end subroutine inter_SmoothVelo_2D

