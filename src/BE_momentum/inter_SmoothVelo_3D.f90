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
! Program unit: inter_SmoothVelo_3D 
! Description: To calculate a corrective term for velocity.    
!-------------------------------------------------------------------------------
subroutine inter_SmoothVelo_3D
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
integer(4) :: npi,i,j,npj,contj,npartint,ibdt,ibdp,ii,Ncbf,icbf,iface,facestr
double precision :: rhoi,rhoj,amassj,pesoj,moddervel,IntWdV,unity
double precision,dimension(1:SPACEDIM) :: DVLoc,DVGlo,BCLoc,BCGlo,LocX
double precision,dimension(3) :: dervel 
double precision,dimension(:),allocatable :: unity_vec
double precision,dimension(:,:),allocatable :: dervel_mat
character(4)     :: strtype
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
   unity = 0.
endif
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
!$omp private(pesoj,Ncbf,icbf,ibdt,ibdp,LocX,iface,facestr,strtype,DVGlo,DVLoc)&
!$omp private(i,j,IntWdV,BCLoc,BCGlo,unity)                                    &
!$omp shared(nag,pg,Med,nPartIntorno,NMAXPARTJ,PartIntorno,PartKernel)         &
!$omp shared(BoundaryDataPointer,BoundaryDataTab,BoundaryFace,Tratto)          &
!$omp shared(indarrayFlu,Array_Flu,esplosione,Domain,n_bodies,unity_vec)       &
!$omp shared(dervel_mat,Granular_flows_options)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   pg(npi)%var = zero
! The mixture particles, which are temporarily affected by the frictional 
! viscosity threshold are fixed.
   if (pg(npi)%mu==Med(pg(npi)%imed)%mumx) cycle
   do contj=1,nPartIntorno(npi)
      npartint = (npi - 1) * NMAXPARTJ + contj
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
      pg(npi)%var(:) = pg(npi)%var(:) + dervel(:) * pesoj   
   enddo
   if (esplosione) pg(npi)%Envar = pg(npi)%Envar + (pg(npj)%IntEn -            &
                                   pg(npi)%IntEn) * pesoj
   if (n_bodies>0) then
      pg(npi)%var(:) = pg(npi)%var(:) + dervel_mat(npi,:)
      unity = unity + unity_vec(npi)
   endif
! Impose boundary conditions at inlet and outlet sections (DB-SPH)
   if (Domain%tipo=="bsph") then
      call DBSPH_inlet_outlet(npi)
      else
         ncbf = BoundaryDataPointer(1,npi)
         ibdt = BoundaryDataPointer(3,npi)
         if (Ncbf>0) then  
            do icbf=1,Ncbf
               ibdp = ibdt + icbf - 1
               LocX(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
               iface = BoundaryDataTab(ibdp)%CloBoNum
               facestr = BoundaryFace(iface)%stretch
               strtype = Tratto(facestr)%tipo
               if ((strtype=='sour').or.(strtype=='velo').or.                  &   
                  (strtype=='flow')) then
                  pg(npi)%var(:) = zero   
                  exit  
               endif
               DVGlo(:) = two * (Tratto(facestr)%velocity(:) - pg(npi)%vel(:))
               do i=1,SPACEDIM
                  DVLoc(i) = zero
                  do j=1,SPACEDIM
                     DVLoc(i) = DVLoc(i) + BoundaryFace(iface)%T(j,i) * DVGlo(j)
                  enddo
               enddo
               IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
               if ((strtype=='fixe').or.(strtype=='tapi')) then
                  BCLoc(1) = DVLoc(1) * IntWdV * Tratto(facestr)%ShearCoeff
                  BCLoc(2) = DVLoc(2) * IntWdV * Tratto(facestr)%ShearCoeff
                  BCLoc(3) = DVLoc(3) * IntWdV
                  do i=1,SPACEDIM
                     BCGlo(i) = zero
                     do j=1,SPACEDIM
                        BCGlo(i) = BCGlo(i) + BoundaryFace(iface)%T(i,j) *     &
                                   BCLoc(j)
                     enddo
                  enddo
                  pg(npi)%var(:) = pg(npi)%var(:) + BCGlo(:)   
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
end subroutine inter_SmoothVelo_3D

