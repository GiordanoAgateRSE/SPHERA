!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: writime2                   
! Description:               
!----------------------------------------------------------------------------------------------------------------------------------

subroutine writime2(ti,tf,nout)
!------------------------
! Modules
!------------------------ 
use I_O_language_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nout
double precision :: ti,tf
logical,save :: first=.true.
double precision :: telat, tcput, tcpup
double precision,save :: tcpus,telas
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
if (first) then
   tcpus = ti
   tcpup = ti
end if
first = .false.
tcput = tf
telat = tf - telas
write(nout,1001) " "
write(nout,1001) cpulbl,tcput,totlbl,tf-ti,przlbl
tcpup = tf
!------------------------
! Deallocations
!------------------------
return
1001 format(1x,a,f20.4,1x,a,f20.4,1x,a)
end subroutine writime2

