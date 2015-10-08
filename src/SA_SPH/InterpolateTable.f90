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
! Program unit: InterpolateTable                                       
! Description: It interpolates values in the array "Table()" with "nrows" rows and "ncols" columns. Independent variables are 
!              in column 0 of Table():
!                 nicols: number of columns of dependent variables to be interpolated
!                 icol(): list of columns of dependent variables to be interpolated
!                 ivalue(): list of the "nicols" interpolated values 
!               (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine InterpolateTable(xval,nicols,icol,ivalue)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nicols,nr0,nc,ic,nr1
double precision :: xval,xval0,csi,deltaX
integer(4),dimension(1:nicols) :: icol
double precision,dimension(1:nicols) :: ivalue
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
deltaX = ktdelta
xval0 = kerneltab(0,0)
nr0 = Int((xval - xval0) / deltaX)
!------------------------
! Statements
!------------------------
if (nr0<=0) then
   do ic=1,nicols
      nc = icol(ic)
      ivalue(ic) = kerneltab(0,nc)
   enddo
   elseif (nr0>=ktrows) then
      do ic=1,nicols
         nc = icol(ic)
         ivalue(ic) = kerneltab(ktrows,nc)
      enddo
      else
         nr1 = nr0 + 1
         csi = (xval - kerneltab(nr0,0)) / deltaX
         do ic=1,nicols
            nc = icol(ic)
            ivalue(ic) = kerneltab(nr0,nc) + csi * (kerneltab(nr1,nc) -        &
               kerneltab(nr0,nc))
         enddo
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine InterpolateTable

