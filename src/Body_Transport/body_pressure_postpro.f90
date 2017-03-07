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
! Program unit: body_pressure_postpro
! Description: Post-processing for body particle pressure.   
!-------------------------------------------------------------------------------
subroutine body_pressure_postpro
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
integer(4) :: npi,j,npartint,npj
double precision :: Sum_W_vol,W_vol,dis,mod_normal,aux
double precision :: dis_vec(3)
integer,dimension(:),allocatable :: wet
double precision,dimension(:),allocatable :: aux_pres
double precision, external :: w
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
allocate(aux_pres(n_body_part))
allocate(wet(n_body_part))
!------------------------
! Initializations
!------------------------
wet = 0
!------------------------
! Statements
!------------------------
! Loop over the body particles (neighbours: body particles at the surface 
! of the same body)
!$omp parallel do default(none)                                                &
!$omp private(Sum_W_vol,npi,npj,dis_vec,dis,mod_normal,W_vol)                  &
!$omp shared(aux_pres,n_body_part,bp_arr,Domain,wet)
do npi=1,n_body_part
   aux_pres(npi) = 0.d0  
   mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
   if (mod_normal>0.d0) then
      Sum_W_vol = 0.d0
! Loop over the body particles at the surface of the same body
      do npj=1,n_body_part
         mod_normal = dsqrt(dot_product(bp_arr(npj)%normal,bp_arr(npj)%normal))
         if ((bp_arr(npi)%body==bp_arr(npj)%body).and.(mod_normal>0.)) then
            dis_vec(:) = bp_arr(npj)%pos(:) - bp_arr(npi)%pos(:) 
            dis = dsqrt(dot_product(dis_vec,dis_vec))
            W_vol = w(dis,Domain%h,Domain%coefke) * bp_arr(npj)%mass
            aux_pres(npi) = aux_pres(npi) + bp_arr(npj)%pres * W_vol
            Sum_W_vol = Sum_W_vol + W_vol 
         endif
      end do
      if (Sum_W_vol>0.000001) aux_pres(npi) = aux_pres(npi) / Sum_W_vol
   endif
enddo
!$omp end parallel do
! Loop over the body particles: update of pressure values 
!$omp parallel do default(none) private(npi) shared(bp_arr,aux_pres,n_body_part)
do npi=1,n_body_part
   if (aux_pres(npi).ne.0.) bp_arr(npi)%pres = aux_pres(npi)
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
deallocate(aux_pres)
deallocate(wet)
return
end subroutine body_pressure_postpro

