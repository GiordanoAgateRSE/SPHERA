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
double precision :: temp_dden,aux,dis,dis_min,x_min,x_max,y_min,y_max,z_min 
double precision :: z_max,mod_normal,W_vol,sum_W_vol
double precision :: dvar(3),aux_vec(3),aux_nor(3),aux_vec2(3)
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
! Loop over the body particles 
! Draft for omp parallelization with critical section 
!omp parallel do default(none)                                                    
!omp private(npi,j,npartint,npj,dvar,temp_dden)                                &
!omp shared(n_body_part,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f,bp_arr)   &
!omp shared(pg,KerDer_bp_f_cub_spl,rag_bp_f)
do npi=1,n_body_part
   bp_arr(npi)%vel_mir = 0.
   sum_W_vol = 0.
! Loop over fluid particles (contributions to fluid particles, similar to a 
! discretized semi-analytic approach or a mirror particle technique) 
   do j=1,nPartIntorno_bp_f(npi)
      npartint = (npi-1)* NMAXPARTJ + j
      npj = PartIntorno_bp_f(npartint)
! Continuity equation
! Normal for u_SA
      dis_min = 999999999.
      x_min = min(bp_arr(npi)%pos(1),pg(npj)%coord(1))
      x_max = max(bp_arr(npi)%pos(1),pg(npj)%coord(1))
      y_min = min(bp_arr(npi)%pos(2),pg(npj)%coord(2))
      y_max = max(bp_arr(npi)%pos(2),pg(npj)%coord(2))
      z_min = min(bp_arr(npi)%pos(3),pg(npj)%coord(3))
      z_max = max(bp_arr(npi)%pos(3),pg(npj)%coord(3))
      aux_nor = 0.
      do k=1,n_body_part
         if (bp_arr(k)%body==bp_arr(npi)%body) then
            aux=dot_product(rag_bp_f(:,npartint),bp_arr(k)%normal)
            if (aux>0.) then
               if (npi==k) then
                  aux_nor(:) = bp_arr(k)%normal(:)
                  exit
                  else
                     if ((bp_arr(k)%pos(1)>=x_min).and.(bp_arr(k)%pos(1)<=     &
                        x_max).and.(bp_arr(k)%pos(2)>=y_min).and.              &
                        (bp_arr(k)%pos(2)<=y_max).and.(bp_arr(k)%pos(3)>=z_min)&
                        .and.(bp_arr(k)%pos(3)<=z_max)) then
                        call distance_point_line_3D                            &
                           (bp_arr(k)%pos,bp_arr(npi)%pos,pg(npj)%coord,dis)
                        dis_min =min(dis_min,dis)
                        if (dis==dis_min) aux_nor(:) = bp_arr(k)%normal(:)
                     endif
               endif
            endif
         endif
      end do
! In case the normal is not yet defined, the normal vector is that of the body
! particle (at the body surface), which is the closest to the fluid particle
      mod_normal = dsqrt(dot_product(aux_nor,aux_nor))
      if (mod_normal==0.) then
         dis_min = 999999999.
         do k=1,n_body_part
            mod_normal=dot_product(bp_arr(k)%normal,bp_arr(k)%normal)
            if (mod_normal>0.) then
               aux_vec(:) = bp_arr(k)%pos(:) - pg(npj)%coord(:) 
               dis = dsqrt(dot_product(aux_vec,aux_vec))
               dis_min = min(dis_min,dis)
               if (dis==dis_min) aux_nor(:) = bp_arr(k)%normal(:)
            endif
         enddo
      endif
! Relative velocity for the continuity equation     
      aux_vec2(:) = bp_arr(npi)%vel(:) - pg(npj)%vel(:)
      dvar(:) = aux_nor(:) * 2. * dot_product(aux_vec2,aux_nor)
      dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
      W_vol = w(dis,Domain%h,Domain%coefke) * (pg(npj)%mass/pg(npj)%dens)
      bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) + (dvar(:)+              &
                               pg(npj)%vel(:)) * W_vol
      sum_W_vol = sum_W_vol + W_vol            
! Contributions to the continuity equation       
      temp_dden = pg(npj)%mass / (dx_dxbodies ** ncord) *                      &
                  KerDer_bp_f_cub_spl(npartint) * (dvar(1) * ( -               &
                  rag_bp_f(1,npartint)) + dvar(2) * ( - rag_bp_f(2,npartint))  &
                  + dvar(3) * ( - rag_bp_f(3,npartint)))
      pg(npj)%dden = pg(npj)%dden - temp_dden
   end do
   if (sum_W_vol>0.) bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) / sum_W_vol     
end do
! Draft for omp parallelization with critical section 
!omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine body_particles_to_continuity

