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
! Program unit: ComputeKernelTable                                   
! Description: To pre-compute and store in kerneltab(0:ktrows,0:ktcols) the following values:
!                 kerneltab(0:ktrows, 0) = rob = rb/h
!                 kerneltab(0:ktrows, 1) = Int W* ro2 dro         (from rob to 2)
!                 kerneltab(0:ktrows, 2) = Int dW*/dro ro dro     (from rob to 2)
!                 kerneltab(0:ktrows, 3) = Int dW*/dro ro^2 dro   (from rob to 2)
!                 kerneltab(0:ktrows, 4) = Int dW*/dro ro^3 dro   (from rob to 2) 
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ComputeKernelTable
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nr
double precision :: rob
character(len=lencard) :: nomsub = "ComputeKernelTable"
double precision,external :: IWro2dro,JdWsRn
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ktdelta = two / INT_KERNELTABLE
!------------------------
! Statements
!------------------------
if (ncord==3) then
   do nr=0,ktrows
     rob = ktdelta * nr
     kerneltab(nr,0) = rob
     kerneltab(nr,1) = IWro2dro(rob)
     kerneltab(nr,2) = JdWsRn(rob,3,1,1) * Unosusquareh
     kerneltab(nr,3) = JdWsRn(rob,3,2,1) * Unosuh
     kerneltab(nr,4) = JdWsRn(rob,3,3,1)
   enddo
   elseif (ncord==2) then
      call diagnostic (arg1=8,arg2=3,arg3=nomsub)
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeKernelTable

