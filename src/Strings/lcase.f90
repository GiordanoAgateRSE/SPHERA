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
! Program unit: lcase                                           
! Description: 
!----------------------------------------------------------------------------------------------------------------------------------
character(100) function lcase(ainp)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: n,ia
character(*) :: ainp
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
do n=1,100
   lcase(n:n) = " "
enddo
!------------------------
! Statements
!------------------------
do n=1,len_trim(ainp)
   ia = iachar(ainp(n:n))
   if ((ia>=65).and.(ia<=90)) then
      lcase(n:n) = char(ia+32)
      else
         lcase(n:n) = ainp(n:n)
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end function lcase

