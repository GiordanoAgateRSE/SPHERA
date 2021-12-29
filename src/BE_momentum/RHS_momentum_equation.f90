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
! Program unit: RHS_momentum_equation
! Description: Right Hand Side of the mometum equation
!-------------------------------------------------------------------------------
subroutine RHS_momentum_equation
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Dynamic_allocation_module
use Hybrid_allocation_module
#ifdef SPACE_3D
use I_O_file_module
use Memory_I_O_interface_module
#endif
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: ii,npi
#ifdef SPACE_3D
integer(4) :: Ncbf_Max,Ncbf
#elif defined SPACE_2D
integer(4) :: Ncbs,IntNcbs
#endif
double precision,dimension(1:SPACEDIM) :: tpres,tdiss,tvisc,BoundReaction
#ifdef SPACE_3D
character(len=lencard) :: nomsub = "RHS_momentum_equation"
#endif
double precision,dimension(1:size(Partz)) :: slip_coeff_counter
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
Partz(:)%avg_comp_slip_coeff = 0.d0
Partz(:)%avg_ni_T_SASPH = 0.d0
Partz(:)%avg_tau_wall_f = 0.d0
#ifdef SPACE_3D
Ncbf_Max = 0
#endif
slip_coeff_counter(:) = 0
!------------------------
! Statements
!------------------------
!$omp parallel do default(none)                                                &
!$omp shared(pg,Domain,BoundaryDataPointer,indarrayFlu,Array_Flu,Med)          &
#ifdef SPACE_3D
!$omp shared(Ncbf_Max)                                                         &
#endif
!$omp shared(nag,slip_coeff_counter)                                           &
!$omp private(npi,ii,tpres,tdiss,tvisc)                                        &
#ifdef SPACE_3D
!$omp private(Ncbf)                                                            &
#elif defined SPACE_2D
!$omp private(Ncbs,IntNcbs)                                                    &
#endif
!$omp private(BoundReaction)
! Loop over particles
do ii = 1,indarrayFlu
   npi = Array_Flu(ii)
! The mixture particles in the elastic-plastic strain regime are held fixed
   if (pg(npi)%mu>(Med(pg(npi)%imed)%mumx*(1.d0-1.d-9))) then
      pg(npi)%acc(:) = zero
      cycle
   endif
! Inner terms of the momentum equation
   call inter_EqMoto(npi,tpres,tdiss,tvisc)
! Searching for the boundary faces/sides, which are the nearest the npi-th 
! current particle
   if ((Domain%time_stage==1).or.(Domain%time_split==1)) then 
      pg(npi)%kodvel = 0
      pg(npi)%velass = zero
   endif
   if (Domain%tipo=="semi") then
#ifdef SPACE_3D
         Ncbf = BoundaryDataPointer(1,npi)
#elif defined SPACE_2D
            Ncbs = BoundaryDataPointer(1,npi)
            IntNcbs = BoundaryDataPointer(2,npi)
#endif
      else
#ifdef SPACE_3D
            Ncbf = 0
#elif defined SPACE_2D
               Ncbs = 0
               IntNcbs = 0
#endif
   endif
   if ((Domain%tipo=="semi").and.                                              &
#ifdef SPACE_3D   
      (Ncbf>0)) then
#elif defined SPACE_2D
      (Ncbs>0).and.(IntNcbs>0)) then
#endif
#ifdef SPACE_3D
!$omp critical (omp_Ncbf_Max)
         Ncbf_Max = max(Ncbf_Max,Ncbf)
!$omp end critical (omp_Ncbf_Max)
         call AddBoundaryContributions_to_ME3D(npi,Ncbf,tpres,tdiss,tvisc,     &
            slip_coeff_counter)
#elif defined SPACE_2D
            call AddBoundaryContributions_to_ME2D(npi,IntNcbs,tpres,tdiss,     &
               tvisc,slip_coeff_counter)
#endif
      if (pg(npi)%kodvel==0) then
#ifdef SPACE_3D
            call AddElasticBoundaryReaction_3D(npi,Ncbf,BoundReaction)
#elif defined SPACE_2D
               call AddElasticBoundaryReaction_2D(npi,IntNcbs,BoundReaction)
#endif
         pg(npi)%acc(:) = tpres(:) + tdiss(:) + tvisc(:) + Domain%grav(:) +    &
                          BoundReaction(:)
         else
            pg(npi)%acc(:) = zero
      endif
      else
         if (Domain%tipo=="semi") then
            pg(npi)%acc(:) = tpres(:) + tdiss(:) + tvisc(:) + Domain%grav(:)
            elseif (Domain%tipo=="bsph") then
               pg(npi)%acc(:) = (tpres(:) + tdiss(:) + tvisc(:)) /             &
                                pg(npi)%Gamma + Domain%grav(:)
         endif
   endif
enddo
!$omp end parallel do
#ifdef SPACE_3D
if (Ncbf_Max>input_any_t%MAXCLOSEBOUNDFACES) then
   write(ulog,"(a,i5,a,i5)") "Increase parameter MAXCLOSEBOUNDFACES from",     &
      input_any_t%MAXCLOSEBOUNDFACES," to ",Ncbf_Max
   call diagnostic(arg1=9,arg2=3,arg3=nomsub)
endif
#endif
! Update of the average slip coefficient for each boundary zone
do ii=1,size(Partz)
   if (slip_coeff_counter(ii)>0) then
      Partz(ii)%avg_comp_slip_coeff = Partz(ii)%avg_comp_slip_coeff /          &
                                      slip_coeff_counter(ii)
      Partz(ii)%avg_ni_T_SASPH = Partz(ii)%avg_ni_T_SASPH /                    &
                                 slip_coeff_counter(ii)
      Partz(ii)%avg_tau_wall_f = Partz(ii)%avg_tau_wall_f /                    &
                                 slip_coeff_counter(ii)
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine RHS_momentum_equation
