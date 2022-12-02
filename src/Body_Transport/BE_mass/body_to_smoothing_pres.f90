!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: body_to_smoothing_pres 
! Description: Contributions of body particles to pressure partial smoothing 
!              and update of the mirror body-particle pressure (Amicarelli et 
!              al., 2015, CAF; Amicarelli et al., 2020, CPC; Amicarelli et al., 
!              2022, IJCFD)    
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
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
integer(4) :: npi,jj,npartint,npj
! Weight for neighbouring body particles ("nb")
double precision :: W_vol_nb
! Weight for neighbouring fluid particles ("nf")
double precision :: W_vol_nf
double precision :: dis,pres_mir,Sum_W_vol_nf,aux_scal
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
!$omp parallel do default(none)                                                &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f)  &
!$omp shared(rag_bp_f,Domain,AppUnity_vec,sompW_vec,body_arr)                  &
!$omp shared(body_minimum_pressure_limiter,body_maximum_pressure_limiter,pg)   &
!$omp private(npi,Sum_W_vol_nf,jj,npartint,npj,dis,W_vol_nb,pres_mir,W_vol_nf) &
!$omp private(aux_scal)
do npi=1,n_body_part
   bp_arr(npi)%pres = 0.d0
   Sum_W_vol_nf = 0.d0
! Loop over the neighbouring fluid particles 
   do jj=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + jj
      npj = PartIntorno_bp_f(npartint)
      dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
      W_vol_nb = w(dis,Domain%h,Domain%coefke) * bp_arr(npi)%volume
!$omp critical (omp_AppUnity_vec)
      AppUnity_vec(npj) = AppUnity_vec(npj) + W_vol_nb
!$omp end critical (omp_AppUnity_vec)
! Mirror body-particle pressure: start
      call body_pressure_mirror_interaction(npi,npj,npartint,pres_mir,W_vol_nf)
! Visibility criterion only for the representative pressure value of the 
! surface body particles (influence on graphics only)
      aux_scal = dot_product(rag_bp_f(:,npartint),bp_arr(npi)%normal)
      if ((aux_scal>1.d-12).or.(.not.(bp_arr(npi)%surface))) then
         bp_arr(npi)%pres = bp_arr(npi)%pres + pres_mir * W_vol_nf
         Sum_W_vol_nf = Sum_W_vol_nf + W_vol_nf
      endif
! Mirror body-particle pressure: end
!$omp critical (omp_sompW_vec)
      sompW_vec(npj) = sompW_vec(npj) + (pres_mir - pg(npj)%pres) * W_vol_nb
!$omp end critical (omp_sompW_vec)
   enddo
! Unique value representative of the body-particle pressure (for body dynamics)  
   if (Sum_W_vol_nf>1.d-3) bp_arr(npi)%pres = bp_arr(npi)%pres / Sum_W_vol_nf
! Body-particle pressure limiters
   if (body_minimum_pressure_limiter) then
      if (bp_arr(npi)%pres<0.d0) bp_arr(npi)%pres = 0.d0 
   endif
   if (body_maximum_pressure_limiter) then
      if (bp_arr(npi)%pres>body_arr(bp_arr(npi)%body)%p_max_limiter) then
         bp_arr(npi)%pres = body_arr(bp_arr(npi)%body)%p_max_limiter
      endif
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine body_to_smoothing_pres
#endif
