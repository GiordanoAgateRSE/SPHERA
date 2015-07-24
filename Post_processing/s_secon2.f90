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
! Program unit: s_secon2                 
! Description:         
!----------------------------------------------------------------------------------------------------------------------------------

subroutine s_secon2(tempo) 
!----------------------------------------------------------------------
!
!     subroutine: S_SECON2 
!
!     DESCRIZIONE:
!       Fornisce il tempo di CPU user+system ed il tempo di ELAPSED.
!
!     AUTORE: F. Riccio, E. Bon
!
!     NOTE:   Si utilizza la funzione di sistema  ETIME    (Alliant)
!             Si utilizza la funzione di sistema  ETIME_   (RISC6000)
!             Si utilizza la funzione di sistema  TIMEF    (RISC6000)
!
!     PARAMETRI DI INPUT:
!
!     PARAMETRI DI OUTPUT:
!       TEMPO(1) : tempo di CPU  all'atto della chiamata (user + system)
!       TEMPO(2) : tempo elapsed all'atto della chiamata.

!
!     VARIABILI UTILIZZATE:
!
!     ROUTINES CHIAMATE:
!
!
!     DATA: 17-11-93
!
!     SUCCESSIVE MODIFICHE:
!
!     DATA    AUTORE      OGGETTO
!     20-12-94 E.Bon      Adattamento a RISC 6000
!
!
!----------------------------------------------------------------------
implicit none

double precision,dimension(2) :: tempo
double precision,dimension(2) :: t2
integer(4)                    :: time

double precision              :: etime

 tempo(1) = etime(t2) 
 tempo(2) = dfloat(time())

return 

end subroutine s_secon2
!---split

