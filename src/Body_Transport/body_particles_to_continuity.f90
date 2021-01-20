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
! Program unit: body_particles_to_continuity
! Description: Contributions of the body particles to the continuity equation.  
!-------------------------------------------------------------------------------
subroutine body_particles_to_continuity
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
integer(4) :: npi,j,npartint,npj,k
double precision :: temp_dden,aux,dis,dis_min,mod_normal,W_vol,sum_W_vol
double precision :: dis_fb_sb,dis_s0_sb,dis_fb_s0
double precision :: dvar(3),aux_vec(3),aux_nor(3)
double precision,external :: w  
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
!$omp shared(KerDer_bp_f_cub_spl,rag_bp_f,pg,Domain,dx_dxbodies)               &
!$omp shared(FSI_free_slip_conditions,ulog,on_going_time_step,n_surf_body_part)&
!$omp shared(surf_body_part)                                                   &
!$omp private(npi,sum_W_vol,W_vol,j,npartint,npj,k,temp_dden,aux,dis,dis_min)  &
!$omp private(mod_normal,dvar,aux_vec,aux_nor,dis_s0_sb,dis_fb_s0,dis_fb_sb)
do npi=1,n_body_part
   bp_arr(npi)%vel_mir(:) = 0.d0
   sum_W_vol = 0.d0
! Loop over fluid particles (contributions to fluid particles) 
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi - 1) * NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
! Continuity equation
! Normal for fluid-body interactions (only for free-slip conditions): start.
      if (FSI_free_slip_conditions.eqv..true.) then
         dis_min = 1.d9
         aux_nor(:) = 0.d0
! Here the array PartIntorno_bp_bp cannot be used as it only refers to 
! neighbouring body particles of other bodies
! Here the array PartIntorno_bp_f cannot be used as usually as it refers to the 
! fluid neighbours of a body particles, not to the solid neighbours of a fluid 
! particle.
         do k=1,n_body_part
            aux_vec(:) = bp_arr(k)%pos(:) - pg(npj)%coord(:) 
            dis_fb_sb = dsqrt(dot_product(aux_vec,aux_vec))
            if (dis_fb_sb<=2.d0*Domain%h) then
! The candidate solid particle "k" (or "sb") has a non-null normal and lies 
! within the kernel support of the fluid particle "npj" (or "fb")
               if (bp_arr(k)%body==bp_arr(npi)%body) then
                  aux=dot_product(rag_bp_f(:,npartint),bp_arr(k)%normal)
                  if (aux>1.d-12) then
! The fluid particle "npj" and the solid particle "k" can "see" each other
                     if (npi==k) then
                        aux_nor(:) = bp_arr(k)%normal(:)
                        exit
                        else
                           aux_vec(:) = bp_arr(npi)%pos(:) - bp_arr(k)%pos(:)
                           dis_s0_sb = dsqrt(dot_product(aux_vec,aux_vec))
                           aux_vec(:) = bp_arr(npi)%pos(:) - pg(npj)%coord(:) 
                           dis_fb_s0 = dsqrt(dot_product(aux_vec,aux_vec))
                           if ((dis_s0_sb<=dis_fb_s0).and.                     &
                              (dis_fb_sb<=dis_fb_s0)) then
! The candidate solid particle "k" (or "sb") lies "between" the fluid particle 
! "npj" (or "fb") and the solid particle "npi" (or "s0"). The candidate lies 
! within the intersection of two equal spheres of radius dis_fb_s0 and centres 
! x_fb and x_s0.
                              call distance_point_line_3D                      &
                               (bp_arr(k)%pos,bp_arr(npi)%pos,pg(npj)%coord,dis)
                              dis_min=min(dis_min,dis)
                              if (dis==dis_min) aux_nor(:) = bp_arr(k)%normal(:)
                           endif
                     endif
                  endif
               endif
            endif
         enddo
         mod_normal = dsqrt(dot_product(aux_nor,aux_nor))
         if (mod_normal<9.99d-1) then
! In case the normal is not yet defined, there is a fluid-solid mass 
! penetration. The normal vector is taken from the neighbouring surface body 
! particle, which is the closest to the neighbouring fluid particle.
            write(ulog,'(3a)') 'The search for the fluid-body interaction ',   &
               'normals detects a fluid-solid mass penetration (in case of ',  &
               'no-slip conditions there is no detection).'
            write(ulog,'(a,i10,a,i10,a,i10,a,i10,a)')                          &
               'Computational body particle: ',npi,                            &
               ' . Neighbouring fluid particle: ',npj,                         &
               ' . Neighbouring body particle: ',k,                            &
               ' . Current time step: ',on_going_time_step,' .'
! Here the array PartIntorno_bp_bp cannot be used as it only refers to 
! neighbouring body particles of other bodies
! Here the array PartIntorno_bp_f cannot be used as it refers to the 
! fluid neighbours of a body particles, not to the solid neighbours of a fluid 
! particle.
            do k=1,n_surf_body_part
               aux_vec(:) = bp_arr(surf_body_part(k))%pos(:) - pg(npj)%coord(:) 
               dis = dsqrt(dot_product(aux_vec,aux_vec))
               dis_min = min(dis_min,dis)
               if (dis==dis_min) aux_nor(:) = bp_arr(k)%normal(:)
            enddo
         endif
      endif
! Normal for fluid-body interactions (only for free-slip conditions): end.
! Relative velocity for the continuity equation     
      aux_vec(:) = bp_arr(npi)%vel(:) - pg(npj)%vel(:)
      if (FSI_free_slip_conditions.eqv..true.) then
! free-slip conditions
         dvar(:) = aux_nor(:) * 2.d0 * dot_product(aux_vec,aux_nor)
         else
! no-slip conditions
            dvar(:) = 2.d0 * aux_vec(:)
      endif
      dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
      W_vol = w(dis,Domain%h,Domain%coefke) * pg(npj)%mass / pg(npj)%dens
      bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) + (dvar(:) +             &
                               pg(npj)%vel(:)) * W_vol
      sum_W_vol = sum_W_vol + W_vol
! Contributions to the continuity equation
      temp_dden = pg(npj)%mass / (dx_dxbodies ** ncord) *                      &
                  KerDer_bp_f_cub_spl(npartint) * (dvar(1) * ( -               &
                  rag_bp_f(1,npartint)) + dvar(2) * ( - rag_bp_f(2,npartint))  &
                  + dvar(3) * ( - rag_bp_f(3,npartint)))
!$omp critical (omp_body_particles_to_continuity)
      pg(npj)%dden = pg(npj)%dden - temp_dden
!$omp end critical (omp_body_particles_to_continuity)
   enddo
   if (sum_W_vol>1.d-18) bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) /     &
                                                  sum_W_vol
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine body_particles_to_continuity
