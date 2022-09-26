!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: proxy_normals_for_body_particles
! Description: Computation of the array «proxy_normal_bp_f» of the indices of 
!              the neighbouring surface body particles providing the normal to 
!              the "body-particle - fluid-particle" interactions (only for 
!              no-slip conditions) and the boundary velocity (for free-slip and 
!              no-slip conditions).  
!-------------------------------------------------------------------------------
#ifdef SOLID_BODIES
subroutine proxy_normals_for_body_particles
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: npi,jj,npartint,npj,kk,npartint_2,npk
double precision :: aux_scal,dis,dis_min,dis_s0_sb,dis_fb_s0
double precision :: aux_vec(3)
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
! Loop over the body particles
!$omp parallel do default(none)                                                &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f)  &
!$omp shared(rag_bp_f,pg,Domain,proxy_normal_bp_f,ulog,uerr,closest_f_sbp)     &
!$omp shared(on_going_time_step,dis_f_sbp,nPartIntorno_f_sbp,PartIntorno_f_sbp)&
!$omp private(npi,jj,npartint,npj,kk,aux_scal,dis,dis_min)                     &
!$omp private(aux_vec,dis_s0_sb,dis_fb_s0,npartint_2,npk)
do npi=1,n_body_part
! Loop over the neighbouring fluid particles 
   do jj=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + jj
      npj = PartIntorno_bp_f(npartint)
      proxy_normal_bp_f(npartint) = 0
      dis_min = 1.d9
! Check if the normal of the interacting body particle is suitable: start
      aux_scal = dot_product(rag_bp_f(:,npartint),bp_arr(npi)%normal)
! Visibility criterion and check on non-null normal
      if (aux_scal>1.d-12) then
! The computational body particle is a surface body particle with a normal 
! suitable for this interaction
         proxy_normal_bp_f(npartint) = npi
! Check if the normal of the interacting body particle is suitable: end
         else
! Search for a proxy normal from a representative surface body particle
            dis_fb_s0 = dsqrt(dot_product(rag_bp_f(:,npartint),                &
                        rag_bp_f(:,npartint)))
! Loop over the neighbouring surface body particles.
! Here the array "PartIntorno_bp_bp" cannot be used as elsewhere as it only 
! refers to neighbouring body particles of other bodies
! Here the array "PartIntorno_bp_f" cannot be used as elsewhere as it refers to 
! the fluid neighbours of a body particles, not to the solid neighbours of a 
! fluid particle.
            do kk=1,nPartIntorno_f_sbp(npj)
               npartint_2 = (npj - 1) * NMAXPARTJ + kk
               npk = PartIntorno_f_sbp(npartint_2)
! The candidate solid particle "npk" (or "sb") lies 
! within the kernel support of the fluid particle "npj" (or "fb").
! Herefater an explicit exclusion of the interacting body particle would slow 
! down the algorithm; it is better to execute few lines just once without 
! effects.
               if (bp_arr(npk)%body==bp_arr(npi)%body) then
! The proxy normal does not belong to another body to be suitable for the 
! visibility criterion
                  aux_scal = dot_product(rag_bp_f(:,npartint),                 &
                             bp_arr(npk)%normal)
                  if (aux_scal>1.d-12) then
! The fluid particle "npj" and the solid particle "npk" can "see" each other
                     aux_vec(1:3) = bp_arr(npi)%pos(1:3) - bp_arr(npk)%pos(1:3)
                     dis_s0_sb = dsqrt(dot_product(aux_vec,aux_vec))
                     if ((dis_s0_sb<=dis_fb_s0).and.                           &
                        (dis_f_sbp(npartint_2)<=dis_fb_s0)) then
! The candidate solid particle "npk" (or "sb") lies "between" 
! the fluid particle "npj" (or "fb") and the solid particle "npi" (or "s0"): 
! the candidate lies within the intersection of two equal spheres of radius 
! "dis_fb_s0" and centres "x_fb" and "x_s0".
                        call distance_point_line_3D(bp_arr(npk)%pos,           &
                           bp_arr(npi)%pos,pg(npj)%coord,dis)
                        if (dis<dis_min) then
                           dis_min = dis
                           proxy_normal_bp_f(npartint) = npk
                        endif
                     endif
                  endif
               endif
            enddo
            if (proxy_normal_bp_f(npartint)==0) then
! In case the normal is not yet defined, this means that there is a fluid-solid 
! mass penetration or some normals are wrongly defined or dx/dx_s is very large 
! so that not all the surface body particles have been considered above. Under 
! these cimrcumstances, the normal vector is taken from the neighbouring 
! surface body particle, which is the closest to the computational body 
! particle (default choice).
               if (closest_f_sbp(npj)>0) then
                  proxy_normal_bp_f(npartint) = closest_f_sbp(npj)
                  else
! In case the normal is not yet defined, this means that there is a fluid-solid 
! mass penetration or some normals are wrongly defined, so that a fluid 
! particle migh have neighbouring body particles, but not neighbouing surface 
! body particles. Under these cimrcumstances, the normal vector is zeroed 
! (secondary default choice) by considering the computational body particle, 
! which here is not on the surface of its body.
                     proxy_normal_bp_f(npartint) = npi
               endif
            endif
      endif
   enddo
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine proxy_normals_for_body_particles
#endif
