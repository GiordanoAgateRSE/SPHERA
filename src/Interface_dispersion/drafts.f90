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
! Program unit: inter_CoefDif       
! Description: To calculate a corrective term for velocity.                 
!-------------------------------------------------------------------------------
subroutine inter_CoefDif (npi)
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
integer(4),intent(IN) :: npi
integer(4) :: npj,contj,npartint
double precision :: unity,rhoj,amassj,pesoj  
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
unity = zero
pg(npi)%veldif(:) = zero
!------------------------
! Statements
!------------------------
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
   rhoj   = pg(npj)%dens
   amassj = pg(npj)%mass
   pesoj = amassj * PartKernel(4,npartint) / rhoj
   unity = unity + pesoj  
   pg(npi)%veldif(:) = pg(npi)%veldif(:) + pg(npj)%vel(:) * pesoj  
enddo
pg(npi)%uni = unity
!------------------------
! Deallocations
!------------------------
return
end subroutine inter_CoefDif

!-------------------------------------------------------------------------------
! Program unit: inter_SmoothVF        
! Description: To calculate a corrective term for the volume fraction.                  
!-------------------------------------------------------------------------------

subroutine inter_SmoothVF (npi,appo1,unity)
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
integer(4),intent(IN) :: npi
double precision,intent(INOUT) :: appo1, unity
integer(4) :: npj,contj,npartint
double precision :: voli,volj,rhoj,amassj,pesoj
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
unity = zero
appo1 = zero
!------------------------
! Statements
!------------------------
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1)* NMAXPARTJ + contj
   npj = PartIntorno(npartint)
   voli = pg(npi)%VolFra
   rhoj = pg(npj)%dens    
   volj = pg(npj)%VolFra    
   amassj = pg(npj)%mass
   if ( pg(npj)%vel_type/="std") cycle
   pesoj = amassj * PartKernel(4,npartint) / rhoj
   unity = unity + pesoj  
   appo1 = appo1 + (volj - voli ) * pesoj  
end do
pg(npi)%uni = unity
!------------------------
! Deallocations
!------------------------
return
end subroutine inter_SmoothVF

!-------------------------------------------------------------------------------
! Program unit: AggDens        
! Description:                   
!-------------------------------------------------------------------------------

subroutine AggDens
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
integer(4) :: npi,ii
double precision :: tirhoc,tirhow,vdens,a1,b1,c1
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
!$omp parallel do default(none)                                                &
!$omp private(npi,vdens,tirhoc,tirhow,a1,b1,c1,ii)                             &  
!$omp shared(nag,Pg,Med,Domain,dt,on_going_time_step,ncord,indarrayFlu,Array_Flu)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
   if (pg(npi)%koddens/=0) cycle
   vdens  = pg(npi)%dden / pg(npi)%dens
   if(on_going_time_step<=1) then
      pg(npi)%tiroc = pg(npi)%rhoc * pg(npi)%VolFra
      tirhoc = pg(npi)%tiroc + dt * (pg(npi)%tiroc * vdens + pg(npi)%rhoc *    &
               pg(npi)%diffu - Med(2)%settlingcoef)
      else    
         tirhoc = pg(npi)%tiroc + dt * (pg(npi)%tiroc * vdens + pg(npi)%rhoc * &
                  pg(npi)%diffu - Med(2)%settlingcoef)
         pg(npi)%tiroc = tirhoc     
   endif
   if (tirhoc>=(pg(npi)%dens/pg(npi)%VolFra)) then
      tirhoc = pg(npi)%dens
      tirhow = zero
      pg(npi)%VolFra = VFmx
      pg(npi)%rhoc = pg(npi)%dens
      pg(npi)%mass = pg(npi)%dens * (Domain%dx ** ncord)
      elseif (tirhoc<=zero) then
         tirhoc = zero
         tirhow = pg(npi)%dens
         pg(npi)%rhow = pg(npi)%dens
         pg(npi)%VolFra = VFmn
         pg(npi)%mass = pg(npi)%dens * (Domain%dx ** ncord)
         else
            tirhow = pg(npi)%dens - tirhoc
            a1 = med(2)%den0 * med(2)%celerita * med(2)%celerita - med(1)%den0 &
                 * med(1)%celerita * med(1)%celerita
            b1 = - (med(2)%celerita * med(2)%celerita) * (tirhoc + med(2)%den0)&
                 - (med(1)%celerita * med(1)%celerita) * (tirhow - med(1)%den0)
            c1 = (med(2)%celerita * med(2)%celerita) * tirhoc
            pg(npi)%VolFra = ( - b1 - Dsqrt( b1 * b1 - 4.d0 * a1 * c1)) / (two &
                             * a1)   
            if (pg(npi)%VolFra >= VFmx) then
               pg(npi)%VolFra  = VFmx
               pg(npi)%rhoc = pg(npi)%dens 
               pg(npi)%rhow = Med(1)%den0 
               pg(npi)%mass = pg(npi)%dens * (Domain%dx**ncord)
               elseif (pg(npi)%VolFra<=VFmn) then
                  pg(npi)%VolFra = VFmn
                  pg(npi)%rhoc = Med(2)%den0 
                  pg(npi)%rhow = pg(npi)%dens
                  pg(npi)%mass = pg(npi)%dens * (Domain%dx ** ncord)
                  else
                     pg(npi)%rhoc = tirhoc / pg(npi)%VolFra  
                     pg(npi)%rhow = tirhow / (one - pg(npi)%VolFra)  
                     pg(npi)%mass = (pg(npi)%VolFra * pg(npi)%rhoc + (one -    &
                                    pg(npi)%VolFra) * pg(npi)%rhow) *          &
                                    (Domain%dx ** ncord)
            endif
   endif
enddo       
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine AggDens

