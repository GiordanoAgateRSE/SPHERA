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
! Program unit: FindLine
! Description: Finds extremes of the rectangular frame which contains the boundary mib.
!----------------------------------------------------------------------------------------------------------------------------------

subroutine FindLine (Xmin,Xmax,Nt)
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
integer(4),intent(IN) :: Nt
double precision,intent(INOUT),dimension(SPACEDIM,NumTratti) :: Xmin,Xmax
integer(4) :: i,iv,nf
character(len=lencard) :: nomsub = "FindLine"
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
! Loop on the segments of the line
do iv=Tratto(Nt)%inivertex,(Tratto(Nt)%inivertex+Tratto(Nt)%numvertices-1)
! It checks for the maximum boundary lines storage
   if (iv>NumBVertices) then
      call diagnostic (arg1=10,arg2=2,arg3=nomsub)
   end if
! It evaluates the minimum and maximum coordinates with respect the current 
! vertex nf in the line nt
   nf = BoundaryVertex(iv)
   do i=1,SPACEDIM
      if (Vertice(i,nf)<Xmin(i,Nt)) Xmin(i,Nt)=Vertice(i,nf)
      if (Vertice(i,nf)>Xmax(i,Nt)) Xmax(i,Nt)=Vertice(i,nf)
   end do
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine FindLine

