!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
!              the "body-particle - fluid-particle" interactions.  
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
integer(4) :: npi,jj,npartint,npj,kk
double precision :: aux_scal,dis,dis_min,dis_fb_sb,dis_s0_sb
double precision :: dis_fb_s0
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
!$omp shared(rag_bp_f,pg,Domain,surf_body_part,proxy_normal_bp_f,ulog)         &
!$omp shared(on_going_time_step,n_surf_body_part)                              &
!$omp private(npi,jj,npartint,npj,kk,aux_scal,dis,dis_min)                     &
!$omp private(aux_vec,dis_s0_sb,dis_fb_s0,dis_fb_sb)
do npi=1,n_body_part
! Loop over fluid particles 
   do jj=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + jj
      npj = PartIntorno_bp_f(npartint)
      dis_min = 1.d9
      proxy_normal_bp_f(npartint) = 0
! Here the array PartIntorno_bp_bp cannot be used as it only refers to 
! neighbouring body particles of other bodies
! Here the array PartIntorno_bp_f cannot be used as usually as it refers to the 
! fluid neighbours of a body particles, not to the solid neighbours of a fluid 
! particle.
      do kk=1,n_body_part
         aux_vec(1:3) = bp_arr(kk)%pos(1:3) - pg(npj)%coord(1:3) 
         dis_fb_sb = dsqrt(dot_product(aux_vec,aux_vec))
         if (dis_fb_sb<=2.d0*Domain%h) then
! The candidate solid particle "kk" (or "sb") has a non-null normal and lies 
! within the kernel support of the fluid particle "npj" (or "fb")
            if (bp_arr(kk)%body==bp_arr(npi)%body) then
! The proxy normal cannot belong to another body to be suitable for the 
! visibility criterion
               aux_scal=dot_product(rag_bp_f(:,npartint),bp_arr(kk)%normal)
               if (aux_scal>1.d-12) then
! The fluid particle "npj" and the solid particle "kk" can "see" each other
                  if (npi==kk) then
! The computational body particle has a normal suitable for this interaction: 
! the visibility check was already passed.
                     proxy_normal_bp_f(npartint) = kk
                     exit
                     else
! Otherwise a proxy normal is needed
                        aux_vec(1:3) = bp_arr(npi)%pos(1:3) -                  &
                                       bp_arr(kk)%pos(1:3)
                        dis_s0_sb = dsqrt(dot_product(aux_vec,aux_vec))
                        aux_vec(1:3) = bp_arr(npi)%pos(1:3) -                  &
                                       pg(npj)%coord(1:3) 
                        dis_fb_s0 = dsqrt(dot_product(aux_vec,aux_vec))
                        if ((dis_s0_sb<=dis_fb_s0).and.                        &
                           (dis_fb_sb<=dis_fb_s0)) then
! The candidate solid particle "kk" (or "sb") lies "between" the fluid particle 
! "npj" (or "fb") and the solid particle "npi" (or "s0"). The candidate lies 
! within the intersection of two equal spheres of radius dis_fb_s0 and centres 
! x_fb and x_s0.
                           call distance_point_line_3D(bp_arr(kk)%pos,         &
                              bp_arr(npi)%pos,pg(npj)%coord,dis)
                           if (dis<dis_min) then
                              dis_min = dis
                              proxy_normal_bp_f(npartint) = kk
                           endif
                        endif
                  endif
               endif
            endif
         endif
      enddo
      if (proxy_normal_bp_f(npartint)==0) then
! In case the normal is not yet defined, there is a fluid-solid mass 
! penetration or some normals are wrongly defined. Under these cimrcumstances, 
! the normal vector is taken from the neighbouring surface body 
! particle, which is the closest to the neighbouring fluid particle.
         write(ulog,'(3a)') 'The search for the fluid-body interaction ',      &
            'normals detects a fluid-solid mass penetration or some normals ', &
            'were wrongly calculated by the input mesh (in case of ',          &
            'no-slip conditions there is no detection).'
         write(ulog,'(a,i10,a,i10,a,i10,a,i10,a)')                             &
            'Computational body particle: ',npi,                               &
            ' . Neighbouring fluid particle: ',npj,                            &
            ' . Current time step: ',on_going_time_step,' .'
! Here the array PartIntorno_bp_bp cannot be used as it only refers to 
! neighbouring body particles of other bodies
! Here the array PartIntorno_bp_f cannot be used as it refers to the 
! fluid neighbours of a body particles, not to the solid neighbours of a fluid 
! particle.
         do kk=1,n_surf_body_part
            aux_vec(1:3) = bp_arr(surf_body_part(kk))%pos(1:3) -               &
                           pg(npj)%coord(1:3) 
            dis = dsqrt(dot_product(aux_vec,aux_vec))
            if (dis<dis_min) then
               dis_min = dis
               proxy_normal_bp_f(npartint) = kk
            endif
         enddo
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
