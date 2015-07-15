!cfile s_secon2.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
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

