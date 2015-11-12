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
! Program unit: body_to_smoothing_pres 
! Description: Contributions of body particles to pressure partial smoothing (Amicarelli et al., 2015, CAF).    
!----------------------------------------------------------------------------------------------------------------------------------

subroutine body_to_smoothing_pres(sompW_vec,AppUnity_vec)
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
double precision,dimension(nag),intent(inout) :: sompW_vec,AppUnity_vec
integer(4) :: npi,j,npartint,npj
double precision :: W_vol,dis
double precision, external :: w
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
! Loop over body particles (neighbours: fluid particles)
do npi=1,n_body_part
! Loop over the neighbouring fluid particles 
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
      dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
      W_vol = w(dis,Domain%h,Domain%coefke) * ((Domain%dd / dx_dxbodies) **    &
              ncord)
      AppUnity_vec(npj) = AppUnity_vec(npj) + W_vol
      sompW_vec(npj) = sompW_vec(npj) + (bp_arr(npi)%pres - pg(npj)%pres) *    &
                       W_vol
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine body_to_smoothing_pres

