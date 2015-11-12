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
! Program unit: s_ctime                
! Description:         
!----------------------------------------------------------------------------------------------------------------------------------

subroutine s_ctime(nout)
!------------------------
! Modules
!------------------------ 
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nout
integer(4),dimension(8) :: dat_array
character(LEN=8) :: dat
character(LEN=10) :: ct
character(LEN=5) :: zone
character(LEN=160) :: date_exec
character(LEN=3),dimension(12) :: mesi
data mesi/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov",   &
   "Dec"/
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
call DATE_AND_TIME(dat,ct,zone,dat_array)
date_exec = mesi(dat_array(2))//" "//dat(7:8)//", "//dat(1:4)//          &
            " at "//ct(1:2)//":"//ct(3:4)//":"//ct(5:10)//" "//zone//" GMT"
write(nout,'(a)') trim(date_exec)
!------------------------
! Deallocations
!------------------------
return
end subroutine s_ctime
