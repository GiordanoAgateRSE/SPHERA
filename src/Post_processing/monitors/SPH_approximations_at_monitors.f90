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
! Program unit: SPH_approximations_at_monitors
! Description:  SPH approximations of density, pressure and velocity vector at 
!               a monitoring element. Density is not obtained via EOS to 
!               treat also dense granular flows which are multi-phase.
!-------------------------------------------------------------------------------
subroutine SPH_approximations_at_monitors(pglocal)
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
! Monitor point position
type (TyCtlPoint),intent(inout) :: pglocal
integer(4) :: nceli,ncel,igridi,kgridi,jgridi,irang,krang,jrang,fw,test,mm,npj
integer(4) :: irestocell
#ifdef SOLID_BODIES
integer(4) :: sb
#endif
! Womegaj: product of the kernel function times the volume of the neighbouring 
! element
! pres_consi: ith-order consistency SPH approximation of pressure
! dens_consi: ith-order consistency SPH approximation of density
double precision :: rijlocal,Womegaj,pres_cons0,dens_cons0
! abs_det_thresh: absolute value of the determinant of the inverse of any 
! renormalization matrix. It is a high-pass filter to apply the matrix 
! inversion.
double precision :: abs_det_thresh
! vel_consi: ith-order consistency SPH approximation of the velocity vector
! x_consi_fp: ith-order consistency SPH approximation of the position vector 
! (contributions from fluid particles)
double precision,dimension(3) :: r_vec,vel_cons0,x_cons0_fp,grad_W
! grad_p_consi: ith-order consistency SPH approximation of the pressure gradient
! grad_1_fp: SPH approximation of the gradient of 1 (contributions from fluid 
! particles)
! grad___raw: SPH (raw) approximation of a gradient
double precision,dimension(3) :: grad_p_cons0,grad_p_cons1,grad_p_raw,grad_1_fp
! grad_rho_consi: ith-order consistency SPH approximation of the density 
! gradient
double precision,dimension(3) :: grad_rho_cons0,grad_rho_cons1,grad_rho_raw
! gradWomegaj: product of the kernel gradient times the volume of the 
! neighbouring element
double precision,dimension(3) :: gradWomegaj,aux_vec_2,aux_vec,aux_vec_3
#ifdef SOLID_BODIES
! grad_1_fp_bp: SPH approximation of the gradient of 1 (contributions from 
! fluid and body particles)
! grad_1_fp_sbp: SPH approximation of the gradient of 1 (contributions from 
! fluid particles and surface body particles)
double precision,dimension(3) :: grad_1_fp_bp,grad_1_fp_sbp
! x_consi_fp_bp: ith-order consistency SPH approximation of the position 
! (contributions from fluid particles and body particles)
! x_consi_fp_sbp: ith-order consistency SPH approximation of the position 
! (contributions from fluid particles and surface body particles)
double precision,dimension(3) :: x_cons0_fp_bp,x_cons0_fp_sbp
#endif
! grad_vel_consi(1:3,1): ith-order consistency SPH approximation of u
! grad_vel_consi(1:3,2): ith-order consistency SPH approximation of v
! grad_vel_consi(1:3,3): ith-order consistency SPH approximation of w
double precision,dimension(3,3) :: grad_vel_cons1,grad_vel_cons0,grad_vel_raw
double precision,dimension(3,3) :: aux_mat
integer(4),external :: CellIndices,CellNumber
double precision,external :: w
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine MatrixProduct(AA,BB,CC,nr,nrc,nc)
      implicit none
      integer(4),intent(in) :: nr,nrc,nc
      double precision,intent(in),dimension(nr,nrc) :: AA
      double precision,intent(in),dimension(nrc,nc) :: BB
      double precision,intent(inout),dimension(nr,nc) :: CC
   end subroutine MatrixProduct
   subroutine grad_W_sub(kernel_ID,r_vec,grad_W)
      implicit none
      integer(4),intent(in) :: kernel_ID
      double precision,dimension(3),intent(in) :: r_vec
      double precision,dimension(3),intent(out) :: grad_W
   end subroutine grad_W_sub
   subroutine Matrix_Inversion_3x3(mat,inv,abs_det_thresh,test)
      implicit none
      double precision,intent(in) :: abs_det_thresh
      double precision,intent(in) :: mat(3,3)
      double precision,intent(inout) :: inv(3,3)
      integer(4),intent(inout) :: test
   end subroutine Matrix_Inversion_3x3
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
nceli = pglocal%cella
if (nceli==0) return
irestocell = CellIndices(nceli,igridi,jgridi,kgridi)
dens_cons0 = 0.d0
grad_rho_raw(1:3) = 0.d0
grad_rho_cons0(1:3) = 0.d0
grad_rho_cons1(1:3) = 0.d0
pres_cons0 = 0.d0
grad_p_raw(1:3) = 0.d0
grad_p_cons0(1:3) = 0.d0
grad_p_cons1(1:3) = 0.d0
vel_cons0(1:3) = 0.d0
grad_vel_raw(1:3,1:3) = 0.d0
grad_vel_cons0(1:3,1:3) = 0.d0
grad_vel_cons1(1:3,1:3) = 0.d0
x_cons0_fp(1:3) = 0.d0
abs_det_thresh = 1.d-9
grad_1_fp(1:3) = 0.d0
#ifdef SOLID_BODIES
x_cons0_fp_bp(1:3) = 0.d0
x_cons0_fp_sbp(1:3) = 0.d0
grad_1_fp_bp(1:3) = 0.d0
grad_1_fp_sbp(1:3) = 0.d0
#endif
!------------------------
! Statements
!------------------------
! Loop over the grid cells
do jrang=jgridi-1,jgridi+1
   do irang=igridi-1,igridi+1 
      do krang=kgridi-1,kgridi+1
         ncel = CellNumber(irang,jrang,krang)
         if (ncel==0) cycle
! Contributions from fluid particles
         if (Icont(ncel+1)>Icont(ncel)) then
! Loop over the cell particles
            do mm=Icont(ncel),Icont(ncel+1)-1
               npj = NPartOrd(mm)
               if (pg(npj)%vel_type=="fix") cycle
               r_vec(1:3) = pg(npj)%coord(1:3) - pglocal%coord(1:3)
               rijlocal = dot_product(r_vec,r_vec)
               if (rijlocal>square_doubleh) cycle
               rijlocal = dsqrt(rijlocal) 
               Womegaj = pg(npj)%mass * w(rijlocal,Domain%h,Domain%coefke) /   &
                         pg(npj)%dens
               pglocal%sigma_fp = pglocal%sigma_fp + Womegaj
#ifdef SOLID_BODIES
               pglocal%sigma_fp_bp = pglocal%sigma_fp_bp + Womegaj
               pglocal%sigma_fp_sbp = pglocal%sigma_fp_sbp + Womegaj
#endif
! Contribution to the monitoring-element pressure (0th-order consistency)
               pres_cons0 = pres_cons0 + pg(npj)%pres * Womegaj
! Contribution to the monitoring-element density (0th-order consistency)
               dens_cons0 = dens_cons0 + pg(npj)%dens * Womegaj
! Contribution to the monitoring-element velocity (0th-order consistency)
               vel_cons0(1:3) = vel_cons0(1:3) + pg(npj)%vel(1:3) * Womegaj
               if (input_any_t%C1_monitors) then
! Additional statements for 1st-order consistency SPH approximations
! Contributions to x (0th-odrer consistency)
                  x_cons0_fp(1:3) = x_cons0_fp(1:3) + pg(npj)%coord(1:3) *     &
                                    Womegaj
#ifdef SOLID_BODIES
                  x_cons0_fp_bp(1:3) = x_cons0_fp_bp(1:3) + pg(npj)%coord(1:3) &
                                       * Womegaj
                  x_cons0_fp_sbp(1:3) = x_cons0_fp_sbp(1:3) +                  &
                                        pg(npj)%coord(1:3) * Womegaj
#endif
! grad_W
                  call grad_W_sub(kernel_ID=1,r_vec=r_vec,grad_W=grad_W)
                  gradWomegaj(1:3) = grad_W(1:3) * pg(npj)%mass / pg(npj)%dens
! Contributions to raw SPH approximations of the gradients of the following 
! functions: 1, density, pressure, velocity
                  grad_1_fp(1:3) = grad_1_fp(1:3) - gradWomegaj(1:3)
#ifdef SOLID_BODIES
                  grad_1_fp_bp(1:3) = grad_1_fp_bp(1:3) - gradWomegaj(1:3)
                  grad_1_fp_sbp(1:3) = grad_1_fp_sbp(1:3) - gradWomegaj(1:3)
#endif
                  grad_rho_raw(1:3) = grad_rho_raw(1:3) - pg(npj)%dens *       &
                                      gradWomegaj(1:3)
                  grad_p_raw(1:3) = grad_p_raw(1:3) - pg(npj)%pres *           &
                                    gradWomegaj(1:3)
                  grad_vel_raw(1:3,1) = grad_vel_raw(1:3,1) - pg(npj)%vel(1) * &
                                        gradWomegaj(1:3)
                  grad_vel_raw(1:3,2) = grad_vel_raw(1:3,2) - pg(npj)%vel(2) * &
                                        gradWomegaj(1:3)
                  grad_vel_raw(1:3,3) = grad_vel_raw(1:3,3) - pg(npj)%vel(3) * &
                                        gradWomegaj(1:3)
! Contributions to the inverse of the renormalization matrices
                  aux_vec(1:3) = (pg(npj)%coord(1) - pglocal%coord(1)) *       &
                                 gradWomegaj(1:3)
                  aux_vec_2(1:3) = (pg(npj)%coord(2) - pglocal%coord(2)) *     &
                                   gradWomegaj(1:3)
                  aux_vec_3(1:3) = (pg(npj)%coord(3) - pglocal%coord(3)) *     &
                                   gradWomegaj(1:3)
                  pglocal%B_ren_fp(1:3,1) = pglocal%B_ren_fp(1:3,1) +          &
                                            aux_vec(1:3)
                  pglocal%B_ren_fp(1:3,2) = pglocal%B_ren_fp(1:3,2) +          &
                                            aux_vec_2(1:3)
                  pglocal%B_ren_fp(1:3,3) = pglocal%B_ren_fp(1:3,3) +          &
                                            aux_vec_3(1:3)
#ifdef SOLID_BODIES
                  pglocal%B_ren_fp_bp(1:3,1) = pglocal%B_ren_fp_bp(1:3,1) +    &
                                               aux_vec(1:3)
                  pglocal%B_ren_fp_bp(1:3,2) = pglocal%B_ren_fp_bp(1:3,2) +    &
                                               aux_vec_2(1:3)
                  pglocal%B_ren_fp_bp(1:3,3) = pglocal%B_ren_fp_bp(1:3,3) +    &
                                               aux_vec_3(1:3)
                  pglocal%B_ren_fp_sbp(1:3,1) = pglocal%B_ren_fp_sbp(1:3,1) +  &
                                                aux_vec(1:3)
                  pglocal%B_ren_fp_sbp(1:3,2) = pglocal%B_ren_fp_sbp(1:3,2) +  &
                                                aux_vec_2(1:3)
                  pglocal%B_ren_fp_sbp(1:3,3) = pglocal%B_ren_fp_sbp(1:3,3) +  &
                                                aux_vec_3(1:3)
#endif
               endif
            enddo
         endif
! Contributions from DBSPH semi-particles (wall elements; excluded for 
! 1st-order consistency approximations)
         if ((Domain%tipo=="bsph").and.(DBSPH%n_w>0).and.                      &
            (.not.input_any_t%C1_monitors)) then
            if (Icont_w(ncel+1)>Icont_w(ncel)) then
! Loop over the neighbouring wall particles in the cell
               do fw=Icont_w(ncel),Icont_w(ncel+1)-1
                  npj = NPartOrd_w(fw)
! Relative positions and distances
                  r_vec(1:3) = pg_w(npj)%coord(1:3) - pglocal%coord(1:3)
                  rijlocal = dot_product(r_vec,r_vec)
! The wall element lies within the "kernel support" of the monitoring element
                  if (rijlocal>square_doubleh) cycle
                  rijlocal = dsqrt(rijlocal)
                  Womegaj = pg_w(npj)%mass * w(rijlocal,Domain%h,Domain%coefke)&
                            / pg_w(npj)%dens
                  pglocal%sigma_fp = pglocal%sigma_fp + Womegaj
                  pres_cons0 = pres_cons0 + pg_w(npj)%pres * Womegaj
                  dens_cons0 = dens_cons0 + pg_w(npj)%dens * Womegaj                  
                  vel_cons0(1:3) = vel_cons0(1:3) + pg_w(npj)%vel(1:3) * Womegaj      
               enddo
            endif
         endif
#ifdef SOLID_BODIES
! Contributions from body particles
         if (Icont_bp(ncel+1)>Icont_bp(ncel)) then
! Loop over the neighbouring solid particles in the cell
            do sb=Icont_bp(ncel),Icont_bp(ncel+1)-1
               npj = NPartOrd_bp(sb)
! Relative positions and distances
               r_vec(1:3) = bp_arr(npj)%pos(1:3) - pglocal%coord(1:3)
               rijlocal = dot_product(r_vec,r_vec)
               if (rijlocal>square_doubleh) cycle
! The body particle lies within the "kernel support" of the monitoring element
               rijlocal = dsqrt(rijlocal)
               Womegaj = w(rijlocal,Domain%h,Domain%coefke) * bp_arr(npj)%volume
! Contribution to the Shepard coefficient
               pglocal%sigma_fp_bp = pglocal%sigma_fp_bp + Womegaj
! Contribution to the 0-th order consistency approximation of pressure
               pres_cons0 = pres_cons0 + bp_arr(npj)%pres * Womegaj
               if (bp_arr(npj)%surface) then
! Only surface body particles are involved in the SPH approximation of velocity
                  vel_cons0(1:3) = vel_cons0(1:3) + bp_arr(npj)%vel(1:3)       &
                                   * Womegaj
! Contribution to the Shepard coefficient
                  pglocal%sigma_fp_sbp = pglocal%sigma_fp_sbp + Womegaj
               endif
               if (input_any_t%C1_monitors) then
! Additional statements for 1st-order consistency SPH approximations
! Contribution to x
                  x_cons0_fp_bp(1:3) = x_cons0_fp_bp(1:3) +                    &
                                       bp_arr(npj)%pos(1:3) * Womegaj
! grad_W
                  call grad_W_sub(kernel_ID=1,r_vec=r_vec,grad_W=grad_W)
! Contributions to raw SPH approximations of the gradients of the following 
! functions: 1, pressure
                  gradWomegaj(1:3) = grad_W(1:3) * bp_arr(npj)%volume
                  grad_1_fp_bp(1:3) = grad_1_fp_bp(1:3) - gradWomegaj(1:3)
                  grad_p_raw(1:3) = grad_p_raw(1:3) - bp_arr(npj)%pres *       &
                                    gradWomegaj(1:3)
! Contributions to the inverse of the renormalization matrices
                  aux_vec(1:3) = (bp_arr(npj)%pos(1) - pglocal%coord(1)) *     &
                                 gradWomegaj(1:3)
                  aux_vec_2(1:3) = (bp_arr(npj)%pos(2) - pglocal%coord(2)) *   &
                                   gradWomegaj(1:3)
                  aux_vec_3(1:3) = (bp_arr(npj)%pos(3) - pglocal%coord(3)) *   &
                                   gradWomegaj(1:3)
                  pglocal%B_ren_fp_bp(1:3,1) = pglocal%B_ren_fp_bp(1:3,1) +    &
                                               aux_vec(1:3)
                  pglocal%B_ren_fp_bp(1:3,2) = pglocal%B_ren_fp_bp(1:3,2) +    &
                                               aux_vec_2(1:3)
                  pglocal%B_ren_fp_bp(1:3,3) = pglocal%B_ren_fp_bp(1:3,3) +    &
                                               aux_vec_3(1:3)
                  if (bp_arr(npj)%surface) then
! Only surface body particles are involved in the SPH approximation of velocity
! Contribution to x
                     x_cons0_fp_sbp(1:3) = x_cons0_fp_sbp(1:3) +               &
                                           bp_arr(npj)%pos(1:3) * Womegaj
! Contributions to raw SPH approximations of the gradients of the following 
! functions: 1, velocity
                     grad_1_fp_sbp(1:3) = grad_1_fp_sbp(1:3) - gradWomegaj(1:3)
                     grad_vel_raw(1:3,1) = grad_vel_raw(1:3,1) -               &
                                           bp_arr(npj)%vel(1) * gradWomegaj(1:3)
                     grad_vel_raw(1:3,2) = grad_vel_raw(1:3,2) -               &
                                           bp_arr(npj)%vel(2) * gradWomegaj(1:3)
                     grad_vel_raw(1:3,3) = grad_vel_raw(1:3,3) -               &
                                           bp_arr(npj)%vel(3) * gradWomegaj(1:3)
! Contributions to the renormalization matrices
                     pglocal%B_ren_fp_sbp(1:3,1) = pglocal%B_ren_fp_sbp(1:3,1) &
                                                   + aux_vec(1:3)
                     pglocal%B_ren_fp_sbp(1:3,2) = pglocal%B_ren_fp_sbp(1:3,2) &
                                                   + aux_vec_2(1:3)
                     pglocal%B_ren_fp_sbp(1:3,3) = pglocal%B_ren_fp_sbp(1:3,3) &
                                                   + aux_vec_3(1:3)
                  endif
               endif
! Fluid density cannot be convoluted from solid body particles
            enddo
         endif
#endif
      enddo
   enddo
enddo
if (pglocal%sigma_fp<1.d-1) then
! In case the dicrete Shepard coefficient is smaller than 0.1 the monitor 
! element assumes null pressure, velocity and density. It is a simple filter.
   pglocal%pres = 0.d0
   pglocal%dens = 0.d0
   pglocal%vel(1:3) = 0.d0
   else
! 0th-order consistency SPH approximations: temporary variables
      dens_cons0 = dens_cons0 / pglocal%sigma_fp
      pres_cons0 = pres_cons0 / pglocal%sigma_fp_bp
      vel_cons0(1:3) = vel_cons0(1:3) / pglocal%sigma_fp_sbp
      if (.not.input_any_t%C1_monitors) then
! 0th-order consistency SPH approximations
! rho
         pglocal%dens = dens_cons0
! p
         pglocal%pres = pres_cons0
! u,v,w
         pglocal%vel(1:3) = vel_cons0(1:3)
         else
! Statements for the 1st-order consistency SPH approximations
! 0th-order consistency SPH approximation of grad_rho
            grad_rho_cons0(1:3) = grad_rho_raw(1:3) - dens_cons0 *             &
                                  grad_1_fp(1:3)
! 0th-order consistency SPH approximation of grad_p
#ifdef SOLID_BODIES
            grad_p_cons0(1:3) = grad_p_raw(1:3) - pres_cons0 * grad_1_fp_bp(1:3)
#else
            grad_p_cons0(1:3) = grad_p_raw(1:3) - pres_cons0 * grad_1_fp(1:3)
#endif
! 0th-order consistency SPH approximation of grad_u
#ifdef SOLID_BODIES
            grad_vel_cons0(1:3,1) = grad_vel_raw(1:3,1) - vel_cons0(1) *       &
                                    grad_1_fp_sbp(1:3)
#else
            grad_vel_cons0(1:3,1) = grad_vel_raw(1:3,1) - vel_cons0(1) *       &
                                    grad_1_fp(1:3)
#endif
! 0th-order consistency SPH approximation of grad_v
#ifdef SOLID_BODIES
            grad_vel_cons0(1:3,2) = grad_vel_raw(1:3,2) - vel_cons0(2) *       &
                                    grad_1_fp_sbp(1:3)
#else
            grad_vel_cons0(1:3,2) = grad_vel_raw(1:3,2) - vel_cons0(2) *       &
                                    grad_1_fp(1:3)
#endif
! 0th-order consistency SPH approximation of grad_w
#ifdef SOLID_BODIES
            grad_vel_cons0(1:3,3) = grad_vel_raw(1:3,3) - vel_cons0(3) *       &
                                    grad_1_fp_sbp(1:3)
#else
            grad_vel_cons0(1:3,3) = grad_vel_raw(1:3,3) - vel_cons0(3) *       &
                                    grad_1_fp(1:3)
#endif
! 0th-order consistency SPH approximations of x
            x_cons0_fp(1:3) = x_cons0_fp(1:3) / pglocal%sigma_fp
#ifdef SOLID_BODIES
            x_cons0_fp_bp(1:3) = x_cons0_fp_bp(1:3) / pglocal%sigma_fp_bp
            x_cons0_fp_sbp(1:3) = x_cons0_fp_sbp(1:3) / pglocal%sigma_fp_sbp
#endif
! Matrix inversion to get the renormalization matrices
            call Matrix_Inversion_3x3(pglocal%B_ren_fp,aux_mat,abs_det_thresh, &
               test)
            if (test==1) then
               pglocal%B_ren_fp(1:3,1:3) = aux_mat(1:3,1:3)
               else
                  write(ulog,*) "Monitors. 1st-order consistency replaced by ",&
                     "0th-order consistency. Monitor position: ",              &
                     pglocal%coord(1:3)," . Simulation time: ",simulation_time,&
                     " ."
                  pglocal%B_ren_fp(1:3,1:3) = 0.d0
                  pglocal%B_ren_fp(1,1) = -1.d0
                  pglocal%B_ren_fp(2,2) = -1.d0
                  pglocal%B_ren_fp(3,3) = -1.d0
            endif
#ifdef SOLID_BODIES
            call Matrix_Inversion_3x3(pglocal%B_ren_fp_bp,aux_mat,             &
               abs_det_thresh,test)
            if (test==1) then
               pglocal%B_ren_fp_bp(1:3,1:3) = aux_mat(1:3,1:3)
               else
                  write(ulog,*) "Monitors. 1st-order consistency replaced by ",&
                     "0th-order consistency. Monitor position: ",              &
                     pglocal%coord(1:3)," . Simulation time: ",simulation_time,&
                     " ."
                  pglocal%B_ren_fp_bp(1:3,1:3) = 0.d0
                  pglocal%B_ren_fp_bp(1,1) = -1.d0
                  pglocal%B_ren_fp_bp(2,2) = -1.d0
                  pglocal%B_ren_fp_bp(3,3) = -1.d0
            endif
            call Matrix_Inversion_3x3(pglocal%B_ren_fp_sbp,aux_mat,            &
               abs_det_thresh,test)
            if (test==1) then
               pglocal%B_ren_fp_sbp(1:3,1:3) = aux_mat(1:3,1:3)
               else
                  write(ulog,*) "Monitors. 1st-order consistency replaced by ",&
                     "0th-order consistency. Monitor position: ",              &
                     pglocal%coord(1:3)," . Simulation time: ",simulation_time,&
                     " ."
                  pglocal%B_ren_fp_sbp(1:3,1:3) = 0.d0
                  pglocal%B_ren_fp_sbp(1,1) = -1.d0
                  pglocal%B_ren_fp_sbp(2,2) = -1.d0
                  pglocal%B_ren_fp_sbp(3,3) = -1.d0
            endif
#endif
! grad_rho
            call MatrixProduct(pglocal%B_ren_fp,grad_rho_cons0,                &
               grad_rho_cons1,nr=3,nrc=3,nc=1)
            grad_rho_cons1(1:3) = -grad_rho_cons1(1:3)
! grad_p
#ifdef SOLID_BODIES
            call MatrixProduct(pglocal%B_ren_fp_bp,grad_p_cons0,grad_p_cons1,  &
               nr=3,nrc=3,nc=1)
#else
            call MatrixProduct(pglocal%B_ren_fp,grad_p_cons0,grad_p_cons1,     &
               nr=3,nrc=3,nc=1)
#endif
            grad_p_cons1(1:3) = -grad_p_cons1(1:3)
! grad_u
#ifdef SOLID_BODIES
            call MatrixProduct(pglocal%B_ren_fp_sbp,BB=grad_vel_cons0(1:3,1),  &
               CC=grad_vel_cons1(1:3,1),nr=3,nrc=3,nc=1)
#else
            call MatrixProduct(pglocal%B_ren_fp,BB=grad_vel_cons0(1:3,1),      &
               CC=grad_vel_cons1(1:3,1),nr=3,nrc=3,nc=1)
#endif
            grad_vel_cons1(1:3,1) = -grad_vel_cons1(1:3,1)
! grad_v
#ifdef SOLID_BODIES
            call MatrixProduct(pglocal%B_ren_fp_sbp,BB=grad_vel_cons0(1:3,2),  &
               CC=grad_vel_cons1(1:3,2),nr=3,nrc=3,nc=1)
#else
            call MatrixProduct(pglocal%B_ren_fp,BB=grad_vel_cons0(1:3,2),      &
               CC=grad_vel_cons1(1:3,2),nr=3,nrc=3,nc=1)
#endif
            grad_vel_cons1(1:3,2) = -grad_vel_cons1(1:3,2)
! grad_w
#ifdef SOLID_BODIES
            call MatrixProduct(pglocal%B_ren_fp_sbp,BB=grad_vel_cons0(1:3,3),  &
               CC=grad_vel_cons1(1:3,3),nr=3,nrc=3,nc=1)
#else
            call MatrixProduct(pglocal%B_ren_fp,BB=grad_vel_cons0(1:3,3),      &
               CC=grad_vel_cons1(1:3,3),nr=3,nrc=3,nc=1)
#endif
            grad_vel_cons1(1:3,3) = -grad_vel_cons1(1:3,3)
! delta_x for density approximation
            aux_vec(1:3) = pglocal%coord(1:3) - x_cons0_fp(1:3)
! rho
            pglocal%dens = dens_cons0 + dot_product(grad_rho_cons1,aux_vec)
! delta_x for pressure approximation
#ifdef SOLID_BODIES
            aux_vec(1:3) = pglocal%coord(1:3) - x_cons0_fp_bp(1:3)
#else
            aux_vec(1:3) = pglocal%coord(1:3) - x_cons0_fp(1:3)
#endif
! p
            pglocal%pres = pres_cons0 + dot_product(grad_p_cons1,aux_vec)
! delta_x for velocity approximation
#ifdef SOLID_BODIES
            aux_vec(1:3) = pglocal%coord(1:3) - x_cons0_fp_sbp(1:3)
#else
            aux_vec(1:3) = pglocal%coord(1:3) - x_cons0_fp(1:3)
#endif
! u
            aux_vec_2(1:3) = grad_vel_cons1(1:3,1)
            pglocal%vel(1) = vel_cons0(1) + dot_product(aux_vec_2,aux_vec)
! v
            aux_vec_2(1:3) = grad_vel_cons1(1:3,2)
            pglocal%vel(2) = vel_cons0(2) + dot_product(aux_vec_2,aux_vec)
! w
            aux_vec_2(1:3) = grad_vel_cons1(1:3,3)
            pglocal%vel(3) = vel_cons0(3) + dot_product(aux_vec_2,aux_vec)
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine SPH_approximations_at_monitors
