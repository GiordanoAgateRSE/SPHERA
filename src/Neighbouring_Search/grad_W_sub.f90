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
! Program unit: grad_W_sub
! Description: Gradient of the kernel function. kernel_ID:
!              =1 cubic beta-spline cubic kernel as defined in the paper of 
!                 Monaghan & Lattanzio (1985)
!              =2 anti-cluster cubic kernel as defined in the paper of Gallati 
!                 & Braschi (2000, "L'Acqua")
!-------------------------------------------------------------------------------
subroutine grad_W_sub(kernel_ID,r_vec,grad_W)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: kernel_ID
double precision,dimension(3),intent(in) :: r_vec
double precision,dimension(3),intent(out) :: grad_W
double precision :: dW_by_dq,q_scal,abs_r_vec
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
grad_w(1:3) = 0.d0
!------------------------
! Statements
!------------------------
abs_r_vec = dsqrt(dot_product(r_vec,r_vec))
if (abs_r_vec<(Domain%h*1.d-6)) return
q_scal = abs_r_vec / Domain%h
select case(kernel_ID)
   case(1)
! beta-spline cubic kernel
      if (q_scal<=1.d0) then
         dW_by_dq = 1.d0 / (PIGRECO * (Domain%h ** 3)) * (9.d0 / 4.d0 *        &
                    (q_scal ** 2) - 3.d0 * q_scal)
         elseif (q_scal>=2.d0) then
            dW_by_dq = 0.d0
            else
               dW_by_dq = -3.d0 / (4.d0 * PIGRECO * (Domain%h ** 3)) * ((2.d0 -&
                          q_scal) ** 2)
      endif
   case(2)
! anti-cluster cubic kernel
      if (q_scal<=2.d0) then
         dW_by_dq = -3.d0 * ((2.d0 - q_scal) ** 2) * (15.d0 / (64.d0 * PIGRECO &
                    * (Domain%h ** 3)))
         else
            dW_by_dq = 0.d0
      endif
   case default
endselect
grad_w(1:3) = dW_by_dq / Domain%h * r_vec(1:3) / abs_r_vec
#ifdef SPACE_2D
select case(kernel_ID)
   case(1)
! beta-spline cubic kernel
      grad_w(1:3) = grad_w(1:3) * Domain%h * 10.d0 / 7.d0
   case(2)
! anti-cluster cubic kernel
      grad_w(1:3) = grad_w(1:3) * Domain%h * 4.d0 / 3.d0
   case default
endselect
#endif
!------------------------
! Deallocations
!------------------------
return
end subroutine grad_W_sub
