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
! Program unit: I_O_diagnostic_module            
! Description: To provide global interfaces to the subroutine diagnostic.                     
!----------------------------------------------------------------------------------------------------------------------------------

module I_O_diagnostic_module
interface 
   subroutine diagnostic(arg1,arg2,arg3)
      use Static_allocation_module
      integer(4),intent(in) :: arg1
      integer(4),intent(in),optional :: arg2
      character(LEN=lencard),intent(in),optional :: arg3
   end subroutine 
end interface 
end module I_O_diagnostic_module

