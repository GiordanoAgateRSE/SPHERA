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
! Program unit: ReadCheck                    
! Description:                   
!----------------------------------------------------------------------------------------------------------------------------------

logical function ReadCheck(IoErr,Ier,Nrighe,ainp,listadati,ninp,nout)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: IoErr,Ier,Nrighe,ninp,nout
character(*) :: ainp,listadati
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
if (IoErr==0) then
   Ier = 0
   ReadCheck = .TRUE.
   else
      Ier = 4
      ReadCheck = .FALSE.
      write(nout,"(1x,a)")    ">>>>>>>>>>>>>> Warning:"
      write(nout,"(1x,a)")  
      write(nout,"(1x,a,i5)") "Error reading unit:  ",ninp
      write(nout,"(1x,a,a)")  "Expected data:       ",listadati
      write(nout,"(1x,a,i5,a)")                                                &
         "Last input line read:"//trim(ainp)//"(line number:",Nrighe,")"
endif
!------------------------
! Deallocations
!------------------------
return
end function ReadCheck

