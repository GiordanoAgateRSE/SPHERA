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
! Program unit: AddBoundaryContributions_to_ME3D                                
! Description: 3D SASPH boundary terms for the momentum equation of any 
!              computational particle (Di Monaco et al., 2011, EACFM).
!              SASPH contributions to the renormalization matrix for the 3D 
!              pressure-gradient term (only for the first time step).
!              Inversion of the renormalization matrix in 3D, even in the 
!              absence of SASPH neighbours (only for the first time step).
!              3D SASPH contributions to the ALE1 velocity increment and its 
!              associated term in the momentum equation.
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine AddBoundaryContributions_to_ME3D(npi,Ncbf,tpres,tdiss,tvisc,        &
   slip_coeff_counter)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_file_module
use Neighbouring_Search_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi,Ncbf
double precision,intent(inout),dimension(1:SPACEDIM) :: tpres,tdiss,tvisc
double precision,intent(inout),dimension(1:size(Partz)) :: slip_coeff_counter
integer(4) :: sd,icbf,iface,ibdp,sdj,i,mati,stretch,ix,iy,aux_int,aux_int_2,ii
double precision :: IntdWrm1dV,cinvisci,Monvisc,cinviscmult,dvn
double precision :: FlowRate1,Lb,L,minquotanode,maxquotanode,u_t_0
double precision :: celeri,alfaMon,Mmult,IntGWZrm1dV,slip_coefficient
double precision :: aux_scal
double precision,dimension(1:SPACEDIM) :: vb,vi,dvij,dvij_Mor_SASPH,nnlocal 
double precision,dimension(1:SPACEDIM) :: ViscoMon,ViscoShear,LocPi
double precision,dimension(1:SPACEDIM) :: u_t_0_vector,one_Loc
! Unit vector of the unity vector: direction of (1,1,1)
double precision,dimension(1:SPACEDIM) :: one_vec_dir
! Gravity vector in the local reference system of the SASPH boundary
double precision,dimension(1:SPACEDIM) :: gravity_local
double precision,dimension(1:SPACEDIM) :: B_ren_aux_Loc,B_ren_aux_Glo
double precision,dimension(1:SPACEDIM) :: aux_vec
! Pressure-gradient "grad_p" SASPH term
double precision,dimension(1:SPACEDIM) :: gradpt_SA
! SASPH "grad_p" interaction term (local reference system)
double precision,dimension(1:SPACEDIM) :: gradpt_SA_0w_loc
! SASPH "grad_p" interaction term (global reference system)
double precision,dimension(1:SPACEDIM) :: gradpt_SA_0w
! ALE SASPH term
double precision,dimension(1:SPACEDIM) :: ALEt_SA
! ALE SASPH interaction term (local reference system)
double precision,dimension(1:SPACEDIM) :: ALEt_SA_0w_loc
! ALE SASPH interaction term (global reference system)
double precision,dimension(1:SPACEDIM) :: ALEt_SA_0w
character(4) :: stretchtype
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine wall_function_for_SASPH(u_t_0,d_50,r_0w,slip_coefficient_0w,     &
      ni_T_0w)
      implicit none
         double precision,intent(in) :: u_t_0,d_50,r_0w
         double precision,intent(out) :: slip_coefficient_0w,ni_T_0w
   end subroutine wall_function_for_SASPH
   subroutine MatrixProduct(AA,BB,CC,nr,nrc,nc)
      implicit none
      integer(4),intent(in) :: nr,nrc,nc
      double precision,intent(in),dimension(nr,nrc) :: AA
      double precision,intent(in),dimension(nrc,nc) :: BB
      double precision,intent(inout),dimension(nr,nc) :: CC
   end subroutine MatrixProduct
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
mati = pg(npi)%imed
vi(1:3) = pg(npi)%var(1:3)
ALEt_SA(1:3) = 0.d0
gradpt_SA(1:3) = 0.d0
ViscoMon(:) = zero
ViscoShear(:) = zero
if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
   pg(npi)%kodvel = 0
   pg(npi)%velass(:) = zero
endif
one_vec_dir(1:3) = 1.d0 / dsqrt(3.d0)
!------------------------
! Statements
!------------------------
face_loop: do icbf=1,Ncbf
   ibdp = BoundaryDataPointer(3,npi) + icbf - 1 
   iface = BoundaryDataTab(ibdp)%CloBoNum
   stretch = BoundaryFace(iface)%stretch
   stretchtype = Tratto(stretch)%tipo
   nnlocal(:) = BoundaryFace(iface)%T(:,3)
! Nothing to execute for open boundaries
   if (stretchtype=="open") cycle face_loop
   if (stretchtype=="fixe".or.stretchtype=="tapi") then
! Local coordinates for particle Pi   
      LocPi(:) = BoundaryDataTab(ibdp)%LocXYZ(:) 
      if (LocPi(3)>zero) then                
! The face "iface" interacts with particle Pi
         dvn = zero
         do SD=1,SPACEDIM
            vb(SD) = BoundaryFace(iface)%velocity(SD)
            dvij(SD) = two * (vi(SD) - vb(SD))
            dvn = dvn + BoundaryFace(iface)%T(SD,3) * dvij(SD)
         enddo
! Boundary contribution to the "grad_p term": start
! Local components (first assessment)
         do SD=1,SPACEDIM
            gravity_local(SD) = zero
            do sdj=1,SPACEDIM
               gravity_local(SD) = gravity_local(SD) +                         &
                                   BoundaryFace(iface)%T(sdj,SD) *             &
                                   Domain%grav(sdj)
            enddo
         enddo
! Local components (second assessment)
         call MatrixProduct(BoundaryDataTab(ibdp)%IntGiWrRdV,BB=gravity_local, &
            CC=gradpt_SA_0w_loc,nr=3,nrc=3,nc=1)
! Global components
         call MatrixProduct(BoundaryFace(iface)%T,BB=gradpt_SA_0w_loc,         &
            CC=gradpt_SA_0w,nr=3,nrc=3,nc=1)
! Collecting the contributions to acceleration (added later after the possible 
! inversion of the renormalization matrix)
         gradpt_SA(1:3) = gradpt_SA(1:3) + gradpt_SA_0w(1:3)
! Boundary contribution to the "grad_p term": end
! Boundary contribution to the "ALE term": start
! Local components
         ALEt_SA_0w_loc(1:3) = 2.d0 * pg(npi)%pres / pg(npi)%dens *            &
                               BoundaryDataTab(ibdp)%BoundaryIntegral(4:6)
! Global components
         call MatrixProduct(BoundaryFace(iface)%T,BB=ALEt_SA_0w_loc,           &
            CC=ALEt_SA_0w,nr=3,nrc=3,nc=1)
! Contribution to acceleration
         ALEt_SA(1:3) = ALEt_SA(1:3) + ALEt_SA_0w(1:3)
! Boundary contribution to the "ALE term": end
! Contributions of the neighbouring SASPH frontiers to the inverse of the 
! renormalization matrix for grad_p: start
         if ((on_going_time_step==1).and.(input_any_t%C1_BE)) then
! Local components explicitly depending on the unit vector of the unity vector
            do SD=1,SPACEDIM
               one_Loc(SD) = 0.d0
               do sdj=1,SPACEDIM
                  one_Loc(SD) = one_Loc(SD) + BoundaryFace(iface)%T(sdj,SD) *  &
                     one_vec_dir(sdj)
               enddo
            enddo
! Local components (second assessment)
! "IntGiWrRdV", or equivalently "J_3,w" is always computed using the 
! beta-spline cubic kernel, no matter about the renormalization
            call MatrixProduct(BoundaryDataTab(ibdp)%IntGiWrRdV,BB=one_Loc,    &
               CC=B_ren_aux_Loc,nr=3,nrc=3,nc=1)
! Global components
            call MatrixProduct(BoundaryFace(iface)%T,BB=B_ren_aux_Loc,         &
               CC=B_ren_aux_Glo,nr=3,nrc=3,nc=1)
! Contribution to the renormalization matrix
            do ii=1,3
               pg(npi)%B_ren_gradp(ii,1:3) = pg(npi)%B_ren_gradp(ii,1:3) +     &
                                             B_ren_aux_Glo(ii)
            enddo
         endif
! Contributions of the neighbouring SASPH frontiers to the inverse of the 
! renormalization matrix for grad_p: end
! Boundary contributions to the viscosity terms 
         IntGWZrm1dV = BoundaryDataTab(ibdp)%BoundaryIntegral(7)
         IntdWrm1dV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
         alfaMon = Med(mati)%alfaMon
         if (alfaMon>zero) then
! The volume viscosity term, depending on velocity divergence, can be 
! neglected (otherwise it may cause several problems); 
! further Monaghan's term is activated even for separating particles
            celeri = Med(mati)%celerita
            Monvisc = alfaMon * celeri * Domain%h 
            Mmult = -Monvisc * dvn * IntGWZrm1dV
            ViscoMon(:) = ViscoMon(:) + Mmult * nnlocal(:)
         endif
         if (Partz(Tratto(stretch)%zone)%slip_coefficient_mode==1) then
! Slip coefficient and molecular viscosity from input
            slip_coefficient = Partz(Tratto(stretch)%zone)%BC_shear_stress_input
            cinvisci = pg(npi)%kin_visc
            elseif (Partz(Tratto(stretch)%zone)%slip_coefficient_mode==2) then
! Slip coefficient computed
! Particle tangential (relative) velocity (vector)
! Both "dvn" and "T" are defined with an opposite direction
               u_t_0_vector(:) = 0.5d0 * (dvij(:) - dvn *                      &
                                 BoundaryFace(iface)%T(:,3))
! Particle tangential (relative) velocity (absolute value)
               u_t_0 = dsqrt(dot_product(u_t_0_vector(:),u_t_0_vector(:)))
! To assess the slip coefficient and the turbulent viscosity
               if (CLC_flag.eqv..true.) then
                  aux_int_2 = CellIndices(pg(npi)%cella,ix,iy,aux_int)
                  aux_scal = CLC%z0(ix,iy) * 10.d0
                  call wall_function_for_SASPH(u_t_0,aux_scal,LocPi(3),        &
                     slip_coefficient,cinvisci)
                  else
                     aux_scal =                                                &
                        Partz(Tratto(stretch)%zone)%BC_shear_stress_input *    &
                        10.d0 
                     call wall_function_for_SASPH(u_t_0,aux_scal,LocPi(3),     &
                        slip_coefficient,cinvisci)
               endif
               if (slip_coefficient>1.d-12) then
!$omp critical (avg_slip_coefficient_3D)
! Update of the incremental sum for the slip coefficient
                  Partz(Tratto(stretch)%zone)%avg_comp_slip_coeff =            &
                     Partz(Tratto(stretch)%zone)%avg_comp_slip_coeff +         &
                     slip_coefficient
! Update of the incremental sum for the turbulent viscosity
                  Partz(Tratto(stretch)%zone)%avg_ni_T_SASPH =                 &
                     Partz(Tratto(stretch)%zone)%avg_ni_T_SASPH + cinvisci
! Update of the incremental sum for the wall-function shear stress
                  Partz(Tratto(stretch)%zone)%avg_tau_wall_f =                 &
                     Partz(Tratto(stretch)%zone)%avg_tau_wall_f +              &
                     slip_coefficient * pg(npi)%dens * cinvisci * u_t_0 /      &
                     LocPi(3)
! Update the counter for both the slip coefficient, the turbulent viscosity and 
! the wall-function shear stress
                  slip_coeff_counter(Tratto(stretch)%zone) =                   &
                     slip_coeff_counter(Tratto(stretch)%zone) + 1
!$omp end critical (avg_slip_coefficient_3D)
               endif
         endif
         if (cinvisci>zero) then
            if ((pg(npi)%laminar_flag==1).or.                                  &
               (Tratto(stretch)%laminar_no_slip_check.eqv..false.)) then
! The factor 2 is already present in "dvij"
               cinviscmult = cinvisci * IntdWrm1dV * slip_coefficient
               if ((input_any_t%ALE3).and.(.not.(pg(npi)%p0_neg_ALE))) then
                  dvij_Mor_SASPH(1:3) = 2.d0 * (pg(npi)%vel_fluid(1:3) -       &
                                        BoundaryFace(iface)%velocity(1:3))
                  else
                     dvij_Mor_SASPH(1:3) = 2.d0 * (pg(npi)%var(1:3) -          &
                                           BoundaryFace(iface)%velocity(1:3))
               endif
               ViscoShear(1:3) = ViscoShear(1:3) + cinviscmult *               &
                                 dvij_Mor_SASPH(1:3)
            endif
         endif
      endif
      elseif (stretchtype=="velo") then 
         if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
            pg(npi)%kodvel = 2
            pg(npi)%velass(:) = Tratto(stretch)%NormVelocity * nnlocal(:) 
         endif
         return
         elseif (stretchtype=="flow") then 
            if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
               pg(npi)%kodvel = 2
! Normal component of velocity 
               if (BoundaryFace(iface)%CloseParticles>0) then
                  minquotanode = 9999.d0
                  maxquotanode = const_m_9999
                  do i=1,BoundaryFace(iface)%nodes
                     if (minquotanode>                                         &
                        vertice(3,BoundaryFace(iface)%node(i)%name))           &
                        minquotanode =                                         &
                        vertice(3,BoundaryFace(iface)%node(i)%name)
                     if (maxquotanode<                                         &
                        vertice(3,BoundaryFace(iface)%node(i)%name))           &
                        maxquotanode =                                         &
                        vertice(3,BoundaryFace(iface)%node(i)%name)
                  enddo
                  Lb = BoundaryFace(iface)%CloseParticles_maxQuota -           &
                     minquotanode
                  L = maxquotanode - minquotanode
                  FlowRate1 = Tratto(stretch)%FlowRate * Lb / L
                  Tratto(stretch)%NormVelocity = doubleh * FlowRate1 /         &
                     (BoundaryFace(iface)%CloseParticles * Domain%PVolume)
                  else
                     Tratto(stretch)%NormVelocity = zero
               endif
               pg(npi)%velass(:) = Tratto(stretch)%NormVelocity * nnlocal(:) 
            endif
            return
! Inlet sections
            elseif (stretchtype=="sour") then
               if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
                  pg(npi)%kodvel = 2
                  pg(npi)%velass(:) = Tratto(stretch)%NormVelocity * nnlocal(:) 
               endif
               return
   endif
enddo face_loop
! The renormalization matrix for the pressure-gradient term is inverted 
! just after all its components are collected and just before the 
! 1st-order consistency scheme applies to the summation of all the 
! particle-boundary contributions.
if ((on_going_time_step==1).and.(input_any_t%C1_BE)) then
! Inversion of the renormalization matrix
   call B_ren_gradp_inversion(npi)
endif
if (Ncbf==0) return
! grad_p (renormalization at boundaries): start
if ((input_any_t%C1_BE).and.(.not.(pg(npi)%p0_neg_ALE))) then
! Renormalization of the contributions to the pressure-gradient term
   call MatrixProduct(pg(npi)%B_ren_gradp,BB=gradpt_SA,CC=aux_vec,nr=3,        &
      nrc=3,nc=1)
! Pressure gradient SASPH term: contribution to acceleration (renormalization)
   gradpt_SA(1:3) = -aux_vec(1:3)
endif
! grad_p (renormalization at boundaries): end
! Adding boundary contributions to the momentum equation
! grad_p term and ALE term
tpres(1:3) = tpres(1:3) + ALEt_SA(1:3) + gradpt_SA(1:3)
! Contribution to the ALE velocity increment (here it is still an acceleration)
pg(npi)%dvel_ALE1(1:3) = pg(npi)%dvel_ALE1(1:3) + ALEt_SA(1:3)
! Sub-grid term
tdiss(1:3) = tdiss(1:3) - ViscoMon(1:3)
! For the sign of "ViscoShear", refer to the mathematical model
tvisc(1:3) = tvisc(1:3) + ViscoShear(1:3)
!------------------------
! Deallocations
!------------------------
return
end subroutine AddBoundaryContributions_to_ME3D
#endif
