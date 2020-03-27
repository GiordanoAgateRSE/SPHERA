!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
use I_O_file_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: Ncbf_Max,ii,npi,Ncbf,Ncbs,IntNcbs
double precision,dimension(1:SPACEDIM) :: tpres,tdiss,tvisc,BoundReaction
character(len=lencard) :: nomsub = "RHS_momentum_equation"
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
Ncbf_Max = 0
!$omp parallel do default(none)                                                &
!$omp shared(pg,Domain,BoundaryDataPointer,Ncbf_Max,indarrayFlu,Array_Flu,Med) &
!$omp shared(nag,ncord)                                                        &
!$omp private(npi,ii,tpres,tdiss,tvisc,Ncbf,BoundReaction,Ncbs,IntNcbs)
! Loop over particles
do ii = 1,indarrayFlu
   npi = Array_Flu(ii)
! The mixture particles in the elastic-plastic strain regime are held fixed
   if (pg(npi)%mu>(Med(pg(npi)%imed)%mumx-1.d-12)) then
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
      if (ncord==3) then
         Ncbf = BoundaryDataPointer(1,npi)
         else
            Ncbs = BoundaryDataPointer(1,npi)
            IntNcbs = BoundaryDataPointer(2,npi)
      endif
      else
         if (ncord==3) then
            Ncbf = 0
            else
               Ncbs = 0
               IntNcbs = 0
         endif
   endif
   if ((Domain%tipo=="semi").and.((Ncbf>0).or.((Ncbs>0).and.(IntNcbs>0)))) then
      if (ncord==3) then
!$omp critical (omp_Ncbf_Max)
         Ncbf_Max = max(Ncbf_Max,Ncbf)
!$omp end critical (omp_Ncbf_Max)
         call AddBoundaryContributions_to_ME3D(npi,Ncbf,tpres,tdiss,tvisc)
         else
            call AddBoundaryContributions_to_ME2D(npi,IntNcbs,tpres,tdiss,tvisc)
      endif
      if (pg(npi)%kodvel==0) then
         if (ncord==3) then
            call AddElasticBoundaryReaction_3D(npi,Ncbf,BoundReaction)
            else
               call AddElasticBoundaryReaction_2D(npi,IntNcbs,BoundReaction)
         endif
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
if (Ncbf_Max>Domain%MAXCLOSEBOUNDFACES) then
   write(ulog,"(a,i5,a,i5)") "Increase parameter MAXCLOSEBOUNDFACES from",     &
      Domain%MAXCLOSEBOUNDFACES," to ",Ncbf_Max
   call diagnostic(arg1=9,arg2=3,arg3=nomsub)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine RHS_momentum_equation
