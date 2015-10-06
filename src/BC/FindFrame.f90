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
! Program unit: FindFrame
! Description: It finds extremes of the rectangular frame which contains the boundary mib. 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine FindFrame(Xmin,Xmax,Nt)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: Nt
double precision,intent(INOUT),dimension(SPACEDIM,NumFacce) :: Xmin,Xmax
integer(4) :: i,n,iv,nf,nod
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
do iv=Tratto(Nt)%iniface,(Tratto(Nt)%iniface+Tratto(Nt)%numvertices-1)
   nf = BFaceList(iv)
   do n=1,4
      nod = BoundaryFace(nf)%Node(n)%name
      if (nod<=0) cycle
      do i=1,Ncord
         if (Vertice(i,nod)<Xmin(i,Nt)) Xmin(i,Nt) = Vertice(i,nod)
         if (Vertice(i,nod)>Xmax(i,Nt)) Xmax(i,Nt) = Vertice(i,nod)
      end do
   end do
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine FindFrame

