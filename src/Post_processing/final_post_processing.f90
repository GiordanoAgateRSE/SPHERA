!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: final_log
! Description: Closure of the log file
!-------------------------------------------------------------------------------
subroutine final_log
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Dynamic_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: SpCountot,i_medium,OpCountot,EpCountot,EpOrdGridtot
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
SpCount = 0
OpCount = 0
EpCount = 0
!------------------------
! Statements
!------------------------
if (ulog>0) then
   write(ulog,*) " "
   write(ulog,'(a)')                                                           &
"----------------------------------------------------------------------------------------"
   write(ulog,*) " "
   SpCountot = 0
   do i_medium=1,NMedium
      SpCountot = SpCountot + SpCount(i_medium)
      write(ulog,'(a,i15,a,a)')                                                &
         "Number of source particles        :  SpCount = ",SpCount(i_medium),  &
         " medium ",Med(i_medium)%tipo
   enddo
   write(ulog,'(a,i15)') "Number of total source particles  :  SpCountot = ",  &
      SpCountot
   write(ulog,*) " "
   OpCountot = 0
   do i_medium=1,NMedium
      OpCountot = OpCountot + OpCount(i_medium)
      write(ulog,'(a,i15,a,a)')                                                &
         "Number of outgone particles       :  OpCount = ",OpCount(i_medium),  &
         " medium ",Med(i_medium)%tipo
   enddo
   write(ulog,'(a,i15)')                                                       &
      "Number of total outgone particles :  OpCountot = ",OpCountot
   write(ulog,*) " "
   EpCountot = 0
   do i_medium=1,NMedium
      EpCountot = EpCountot + EpCount(i_medium)
      write(ulog,'(a,i15,a,a)')                                                &
         "Number of escaped particles       :  EpCount = ",EpCount(i_medium),  &
         " medium ",Med(i_medium)%tipo
   enddo
   write(ulog,'(a,i15)')                                                       &
      "Number of total escaped particles :  EpCountot = ",EpCountot
   write(ulog,*) " "
   EpOrdGridtot = 0
   do i_medium=1,NMedium
      EpOrdGridtot = EpOrdGridtot + EpOrdGrid(i_medium)
      write(ulog,'(a,i15,a,a)')                                                &
         "Number of escaped particles (OrdGrid1)       :  EpOrdGrid = ",       & 
         EpOrdGrid(i_medium)," medium ",Med(i_medium)%tipo
   enddo
   write(ulog,'(a,i15)')                                                       &
      "Number of total escaped particles (OrdGrid1) :  EpOrdGridtot = ",       &
      EpOrdGridtot
   write(ulog,*) " "
   write(ulog,*) " "
   write(ulog,*) "Final number of particles:       NAG = ",nag
   if (Domain%tipo=="bsph") write(ulog,*)                                      &
      "Final number of wall particles:       DBSPH%n_w = ",DBSPH%n_w
   write(ulog,*) " "
   write(ulog,'(a)')                                                           &
"----------------------------------------------------------------------------------------"
   write(ulog,*) " "
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine final_log
