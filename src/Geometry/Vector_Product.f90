!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: Vector_Product    
! Description: To return in ww the cross product of vectors uu and vv.           
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Vector_Product(uu,VV,ww,SPACEDIM)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: SPACEDIM
double precision,intent(IN),dimension(SPACEDIM) :: uu,VV
double precision,intent(INOUT),dimension(SPACEDIM) :: ww
integer(4) :: i,j,k
integer(4),dimension(3) :: iseg=(/2,3,1/)
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
do i=1,SPACEDIM
   j = iseg(i)
   k = iseg(j)
   ww(i) = uu(j) * VV(k) - uu(k) * VV(j)
end do
!------------------------
! Deallocations
!------------------------
return
end subroutine Vector_Product

