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
! Program unit: CalcPre 
! Description:  Particle pressure estimation.      
!-------------------------------------------------------------------------------
subroutine CalcPre
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
integer(4) :: npi
double precision :: rhorif,c2,crhorif,wrhorif,wc2,cc2 
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
if (diffusione) then
!$omp parallel do default(none) private(npi,crhorif,wrhorif,wc2,cc2)           &
!$omp shared(nag,Pg,Med)
   do npi=1,nag
      if (pg(npi)%koddens/=0) cycle
      if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
      crhorif = Med(2)%den0
      wrhorif = Med(1)%den0
      wc2 = Med(1)%celerita * Med(1)%celerita
      cc2 = Med(2)%celerita * Med(2)%celerita
      if (pg(npi)%imed==1)then
         pg(npi)%pres = wc2 * (pg(npi)%dens - (crhorif * VFmn + wrhorif * (1 - &
                        VFmn))) 
         elseif (pg(npi)%imed==2) then
         pg(npi)%pres = cc2 * (pg(npi)%dens - (crhorif * VFmx + wrhorif * (1 - &
                        VFmx))) 
      endif
   enddo
!$omp end parallel do
   else
! Is this useless?
      if (on_going_time_step>0) then  
!$omp parallel do default(none) private(npi,rhorif,c2)                         & 
!$omp shared(nag,pg,Domain,Med,esplosione) 
! Loop over all the particles
         do npi=1,nag
! It skips the outgone particles
! It skips the particles with velocity type different from "standard"
            if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
            if (esplosione) then
! Modification for specific internal energy
               pg(npi)%pres = (Med(pg(npi)%imed)%gamma - one) * pg(npi)%IntEn  &
                              * pg(npi)%dens
               else
! It evaluats the new pressure particle value 
                  rhorif = Med(pg(npi)%imed)%den0
                  c2 = Med(pg(npi)%imed)%eps / rhorif
                  if ((pg(npi)%dens-rhorif)/=0.d0) pg(npi)%pres = c2 *         &
                     (pg(npi)%dens - rhorif) + Domain%prif
            endif
         enddo
!$omp end parallel do
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine CalcPre

! Draft subroutine
!subroutine contrmach (k,celmax)
!
!use Static_allocation_module
!use Hybrid_allocation_module
!use Dynamic_allocation_module
!
!implicit none
!
!integer(4) :: k
!
!double precision    :: celmax,amaxnmach,amaxnmach2,amachnumb2
!!double precision   :: vmod2,cel2,amachnumb
!
! amaxnmach  = 0.2d0
! cel2       = celmax*celmax
! amaxnmach2 = amaxnmach*amaxnmach
!
! controllo numero mach - scrive su nout
! vmod2 = pg(k)%vel(1)*pg(k)%vel(1) + pg(k)%vel(3)*pg(k)%vel(3)
! amachnumb2 = vmod2 / cel2
!
! if ( amachnumb2 > amaxnmach2 )then    ! -------------------------------
!    amachnumb = Dsqrt(amachnumb2)
!!   pg(k)%vel(:)=(amaxnmach*celmax/Dsqrt(vmod2))*pg(k)%vel(:) 
! end if        ! -------------------------------
!
!return
!end

