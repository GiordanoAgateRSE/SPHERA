!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: JdWsRn                                          
! Description: (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
double precision function JdWsRn(ro,SD,n,kernel)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
double precision,parameter :: KS2D = 0.454728408833987d0 != 10 / (7*pigreco)
double precision,parameter :: KA2D = 0.09947192d0        !=5./(16*pigreco)
double precision,parameter :: KS3D = 0.31830989d0        != 1 / pigreco
double precision,parameter :: KA3D = 0.07460394d0        !=15./(16*pigreco)
integer(4),intent(IN) :: SD, n, kernel
double precision,intent(IN) :: ro
double precision :: ro2,ro3,ro4,ro5,duemro
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ro2 = ro * ro
ro3 = ro2 * ro
ro4 = ro2 * ro2
ro5 = ro4 * ro
JdWsRn = zero
!------------------------
! Statements
!------------------------
select case (kernel)
   case (1)              
! Cubic spline  Kernel
      select Case (n)
         case (0)      
! n = 0
            if ((ro>=zero).and.(ro<one)) then
               JdWsRn = -one + (1.5d0 - 0.75d0 * ro) * ro2
               elseif ((ro>=one).and.(ro<two)) then
                  duemro = two - ro
                  JdWsRn = -0.25d0 * duemro * duemro * duemro
            endif
         case (1)      
! n = 1
            if ((ro>=zero).and.(ro<one)) then
               JdWsRn = -0.75d0 + (one - 0.5625d0 * ro) * ro3
               elseif ((ro>=one).and.(ro<two)) then
                  JdWsRn = -one + (1.5d0 - one * ro + 0.1875d0 * ro2) * ro2
            endif
         case (2)      
! n = 2
            if ((ro>=zero).and.(ro<one)) then
               JdWsRn = -0.7d0 + (0.75d0 - 0.45d0 * ro) * ro4
               elseif (ro>=one .And. ro<two) then
                  JdWsRn = -0.8d0 + (one - 0.75d0 * ro + 0.15d0 * ro2) * ro3
            endif
         case (3)      
! n = 3
            if ((ro>=zero).and.(ro<one)) then
               JdWsRn = -0.75d0 + (0.6d0 - 0.375d0 * ro) * ro5
               elseif ((ro>=one).and.(ro<two)) then
                  JdWsRn = -0.8d0 + (0.75d0 - 0.6d0 * ro + 0.125d0 * ro2) * ro4
            endif
         case default
            JdWsRn = zero
      end Select
      select case (SD)
         case (2)      
! 2D geometry
            JdWsRn = JdWsRn * KS2D
         case (3)      
! 3D geometry
            JdWsRn = JdWsRn * KS3D
         case default
            JdWsRn = zero
      end Select
   case (2)              
! Gallati anti-cluster kernel 
      select case (n)
         case (1)      
! n = 1
            if (ro<two) then
               JdWsRn = -4.0d0 + (6.0d0 - 4.0d0 * ro + 0.75d0 * ro2) * ro2
            endif
         case (2)      
! n = 2
            if (ro<two) then
               JdWsRn = -3.2d0 + (4.0d0 - 3.0d0 * ro + 0.6d0* ro2) * ro3
            endif
         case (3)      
! n = 3
            if (ro<two) then
               JdWsRn = -3.2d0 + (3.0d0 - 2.4d0 * ro + 0.5d0 * ro2) * ro4
            endif
         case default
            JdWsRn = zero
      endselect
      select case (SD)
         case (2)      
! 2D geometry
            JdWsRn = JdWsRn * KA2D
         case (3)      
! 3D geometry
            JdWsRn = JdWsRn * KA3D
         case default
            JdWsRn = zero
      end Select
   case default
end Select
!------------------------
! Deallocations
!------------------------
return
end function JdWsRn

