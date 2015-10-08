!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: time_integration                                          
! Description: Explicit Runge-Kutta time integration schemes. 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine time_integration  
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: ii,npi
double precision :: appo1,appo2,appo3
character(len=lencard) :: nomsub = "time_integration"
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
! Explicit Runge-Kutta time integration schemes (Euler-RK1-, Heun-RK2)
call start_and_stop(2,17)
select case (Domain%RKscheme)
   case (1) 
      call Euler
   case (2) 
      if (Domain%time_stage==1) then
         call Euler
         else
            call Heun
      endif
end select
call start_and_stop(3,17)
! Diffusion coefficient update
if (diffusione) then
   call start_and_stop(2,15)
!$omp parallel do default(none)                                                &
!$omp private(npi,ii,appo1,appo2,appo3)                                        &
!$omp shared(nag,Pg,Med,indarrayFlu,Array_Flu)
   do ii=1,indarrayFlu
      npi = Array_Flu(ii)
      if ((pg(npi)%VolFra==VFmx).and.                                          &
         (pg(npi)%visc==Med(pg(npi)%imed)%mumx/pg(npi)%dens)) then
         pg(npi)%coefdif = zero
         else
            call inter_CoefDif(npi)
            if (pg(npi)%uni>zero) pg(npi)%veldif = pg(npi)%veldif / pg(npi)%uni
            appo1 = (pg(npi)%veldif(1) - pg(npi)%var(1)) * (pg(npi)%veldif(1)  &
               - pg(npi)%var(1))
            appo2 = (pg(npi)%veldif(2) - pg(npi)%var(2)) * (pg(npi)%veldif(2)  &
               - pg(npi)%var(2))
            appo3 = (pg(npi)%veldif(3) - pg(npi)%var(3)) * (pg(npi)%veldif(3)  &
               - pg(npi)%var(3))
            pg(npi)%coefdif = pg(npi)%coefdif * Dsqrt(appo1 + appo2 + appo3)
      endif
   enddo
!$omp end parallel do
   call start_and_stop(3,15)
endif
if (diffusione) then
   call start_and_stop(2,16)
   call aggdens
   call start_and_stop(3,16)
endif
! Equation of State 
call start_and_stop(2,13)
call calcpre  
call start_and_stop(3,13)
!------------------------
! Deallocations
!------------------------
return
end subroutine time_integration

