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
! Program unit: ltrim                                           
! Description: 
!----------------------------------------------------------------------------------------------------------------------------------

character(10) function ltrim(txt)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,l,n
character(*) :: txt
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
l = len_trim(txt)
!------------------------
! Statements
!------------------------
do n=1,l
   if (txt(n:n)/=" ") then
      txt(1:l-n+1) = txt(n:l)
      do i=(l-n+2),l
         txt(i:i) = " "
      enddo
      exit
   endif
enddo
ltrim = trim(txt)
!------------------------
! Deallocations
!------------------------
return
end function ltrim

