!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: AddBoundaryContribution_to_CE3D                                 
! Description: To compute boundary terms for the 3D continuity equation (rodivV). Equation refers to particle npi. 
!              It performs implicit computation of gradPsuro. 
!              (Di Monaco et al., 2011, EACFM).                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine AddBoundaryContribution_to_CE3D(npi,Ncbf,BCtorodivV)
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
integer(4),intent(IN)    :: npi,Ncbf
double precision,intent(INOUT) :: BCtorodivV
integer(4) :: sd,sdj,icbf,iface,ibdt,ibdp,stretch        
double precision :: roi,scaprod
double precision,dimension(1:SPACEDIM) :: vb,vi,dvij,LocPi, LocDvij
character(4) :: boundtype
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
roi = pg(npi)%dens
vi(:) = pg(npi)%var(:)
BCtorodivV = zero
if (Ncbf<=0) return
ibdt = BoundaryDataPointer(3,npi)
!------------------------
! Statements
!------------------------
do icbf=1,Ncbf
   ibdp = ibdt + icbf - 1
   LocPi(1:SPACEDIM) = BoundaryDataTab(ibdp)%LocXYZ(1:SPACEDIM)
   iface = BoundaryDataTab(ibdp)%CloBoNum
   stretch = BoundaryFace(iface)%stretch
   boundtype = Tratto(stretch)%tipo
   if (boundtype=="fixe".OR.boundtype=="tapi") then
! The face "iface" interacts with particle Pi
      if (LocPi(3)>zero) then          
         vb(:) = BoundaryFace(iface)%velocity(:)
         dvij(:) = two * (vi(:) - vb(:))
         scaprod = zero
         sd = 3
! Local components of the vector "2*(vi-vb)"
         Locdvij(sd) = zero  
         do sdj=1,SPACEDIM
            Locdvij(sd) = Locdvij(sd) + dvij(sdj) *                            &
               BoundaryFace(iface)%T(sdj,sd)
         enddo
         scaprod = scaprod + Locdvij(sd) *                                     &
            BoundaryDataTab(ibdp)%BoundaryIntegral(3+sd)        
! Boundary contribution to the continuity equation 
         BCtorodivV = BCtorodivV + roi * scaprod
      endif
      elseif (boundtype=="velo".or.boundtype=="flow".or.boundtype=="sour") then
         if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
            pg(npi)%koddens = 2
            pg(npi)%densass = roi
            BCtorodivV = zero
         endif
         return
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContribution_to_CE3D

