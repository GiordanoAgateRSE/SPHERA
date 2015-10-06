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
! Program unit: GetToken                                           
! Description: 
!----------------------------------------------------------------------------------------------------------------------------------

character(80) function GetToken(ainp,itok,ioerr)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: itok,ioerr
character(*) :: ainp
logical :: blank
integer(4) :: n,number_token
integer(4),dimension(2,20) :: index_token
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
number_token = 0
index_token  = 0
blank = .TRUE.
!------------------------
! Statements
!------------------------
do n=1,len_trim(ainp)
   if ((blank).and.(ainp(n:n)/=" ")) then
      number_token = number_token + 1
      index_token(1,number_token) = n
      index_token(2,number_token) = n
      blank = .FALSE.
      elseif ((.not.blank).and.(ainp(n:n)/=" ")) then
         index_token(2,number_token) = n
         elseif (ainp(n:n)==" ") then 
         blank = .TRUE.
   endif
end do
if (itok<=number_token) then
   ioerr = 0
   GetToken = ainp(index_token(1,itok):index_token(2,itok))
   else
   ioerr = itok
   GetToken = ""
endif
!------------------------
! Deallocations
!------------------------
return
end function GetToken

