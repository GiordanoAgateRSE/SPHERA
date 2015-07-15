!cfile OrdGrid1.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : OrdGrid1
!
! Last updating : April 18, 2013
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Nov11        BSPH: wall element ordering
! 04  Amicarelli/Agate  13nov12        (AA501b) Body dynamics
! 05  Amicarelli/Agate  18apr13        add check on body_part_reorder
!
!************************************************************************************
! Module purpose : Module for sorting on the particles grid
!
! Calling routine: Gest_Input, Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: ncella
!
!************************************************************************************
!
subroutine OrdGrid1 ( nout )
!* ordinamento su griglia delle particelle
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
integer(4),intent(IN) :: nout
!
!.. Local Scalars ..
integer(4) :: npi,ncel,i
!
!.. Local Arrays ..
integer(4),dimension(Grid%nmax) :: numpartincelgiaaposto
!
!.. External Routines ..
integer(4),external :: ParticleCellNumber
!
!.. Executable Statements ..
!
!* annullamento ordine preesistente
 numpartincelgiaaposto = 0
 Icont                 = 0
!
!* 1 passata
!* si trova la cella in cui cade la particella e si contano le
!* particelle che ci sono in ogni cella
!
!!!!!!  do npi = 1,nag
!!!!!!    if ( pg(npi)%cella == 0 ) cycle
!!!!!!    ncel = ParticleCellNumber(npi)
!!!!!!    if ( ncel <= 0 .or. ncel > Grid%nmax) then
!!!!!!!
!!!!!!!.. the particle has been detected out of the grid and it is removed from the particle array
!!!!!!!
!!!!!!      write(nout,'(a,i7,a,i7,3x,3e15.8)') "ORDGRID1 particle #",npi, "   cell:",ncel,pg(npi)%coord(:)
!!!!!!!ag      EpOrdGrid(pg(npi)%imed) = EpOrdGrid(pg(npi)%imed) + 1
!!!!!!      pg(npi)%cella = 0
!!!!!!!ag      pg(npi) = PgZero
!!!!!!    else
!!!!!!!
!!!!!!!.. count the particles in each cell and the total number of particles (location nmax+1)
!!!!!!!
!!!!!!      Icont(ncel) = Icont(ncel) + 1
!!!!!!      pg(npi)%cella = ncel
!!!!!!      Icont(Grid%nmax+1) = Icont(Grid%nmax+1) + 1

!!!!!!    end if
!!!!!!  end do
!
!ag
!.. inserire eliminazione della particella e fare loop con while
!.. togliere EpOrdGrid(pg(npi)%imed) = EpOrdGrid(pg(npi)%imed) + 1    e   pg(npi) = PgZero   sopra
  npi = 0
  do while (nag > 0 .AND. npi < nag)
    npi = npi + 1
    ncel = ParticleCellNumber(pg(npi)%coord)
    if (pg(npi)%cella <= 0) then
      pg(npi) = pg(nag)
      pg(nag) = PgZero
!
!AA402
     if (Domain%RKscheme > 1) ts0_pg(nag) = ts_pgZero
!
      nag = nag - 1
      npi = npi - 1
!
    else if (ncel <= 0 .or. ncel > Grid%nmax) then
!
!.. the particle has been detected out of the grid and it is removed from the particle array
      write(nout,'(a,i7,a,i7,3x,3e15.8,a,i7)') &
            "ORDGRID1 particle #",npi, "   cell:",ncel,pg(npi)%coord(:),"  Total number of particles is = ",nag-1
      EpOrdGrid(pg(npi)%imed) = EpOrdGrid(pg(npi)%imed) + 1
      pg(npi) = pg(nag)
      pg(nag) = PgZero
!
!AA402
      if (Domain%RKscheme > 1) ts0_pg(nag) = ts_pgZero
!
      nag = nag - 1
      npi = npi - 1
!
    else
!
!.. count the particles in each cell and the total number of particles (location nmax+1)
      Icont(ncel) = Icont(ncel) + 1
      pg(npi)%cella = ncel
      Icont(Grid%nmax+1) = Icont(Grid%nmax+1) + 1
    end if
  end do
!
!ag
!
!* 2 passata
!* si definisce il puntatore d'inizio cella
  do i = Grid%nmax, 1, -1
    Icont(i) = Icont(i+1) - Icont(i)
  end do
!
!* 2 passata bis
!* si incrementa il contatore di 1
  do i = 1, Grid%nmax+1
    Icont(i) = Icont(i) + 1
  end do
!
!* 3 passata finale
!* si mettono le particelle al loro posto nel vettore NPartOrd
!
 do npi = 1,nag
    ncel = pg(npi)%cella
    if ( ncel == 0 ) cycle
    NPartOrd ( Icont(ncel) + numpartincelgiaaposto(ncel) ) = npi
    numpartincelgiaaposto(ncel) = numpartincelgiaaposto(ncel) + 1
 end do
 !
!AA406 start
!AA503sub 
  if ((DBSPH%n_w > 0) .and. ((it_corrente == it_start) .or. (Domain%body_part_reorder == 1)) ) then
!
    Icont_w = 0
    numpartincelgiaaposto = 0
!
!* 1 passata
!* si trova la cella in cui cade la particella e si contano le
!* particelle che ci sono in ogni cella
!
    npi = 0
!AA601 sub    
    do while (npi < (DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet))
       npi = npi + 1
       pg_w(npi)%cella = ParticleCellNumber(pg_w(npi)%coord)
!.. count the particles in each cell and the total number of particles (location nmax+1)
       Icont_w(pg_w(npi)%cella) = Icont_w(pg_w(npi)%cella) + 1
       Icont_w(Grid%nmax+1) = Icont_w(Grid%nmax+1) + 1
    end do
!
!* 2 passata
!* si definisce il puntatore d'inizio cella
    do i = Grid%nmax, 1, -1
       Icont_w(i) = Icont_w(i+1) - Icont_w(i)
    end do
!
!* 2 passata bis
!* si incrementa il contatore di 1
    do i = 1, Grid%nmax+1
       Icont_w(i) = Icont_w(i) + 1
    end do
!
!* 3 passata finale
!* si mettono le particelle al loro posto nel vettore NPartOrd
!AA601
    do npi = 1,(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)
       NPartOrd_w(Icont_w(pg_w(npi)%cella)+numpartincelgiaaposto(pg_w(npi)%cella)) = npi
       numpartincelgiaaposto(pg_w(npi)%cella) = numpartincelgiaaposto(pg_w(npi)%cella) + 1
!AA501
!       pg_w(npi)%wet = 1
!
    end do
! 
 endif
!AA406 end
!

!AA501b start
 if (n_bodies > 0) then
    Icont_bp = 0
    numpartincelgiaaposto = 0

!* 1 passata
!* si trova la cella in cui cade la particella e si contano le
!* particelle che ci sono in ogni cella
!
    npi = 0
    do while (npi < n_body_part)
       npi = npi + 1
       bp_arr(npi)%cell = ParticleCellNumber(bp_arr(npi)%pos)
!.. count the particles in each cell and the total number of particles (location nmax+1)
       Icont_bp(bp_arr(npi)%cell) = Icont_bp(bp_arr(npi)%cell) + 1
       Icont_bp(Grid%nmax+1) = Icont_bp(Grid%nmax+1) + 1
    end do

!* 2 passata
!* si definisce il puntatore d'inizio cella
    do i = Grid%nmax, 1, -1
       Icont_bp(i) = Icont_bp(i+1) - Icont_bp(i)
    end do
!
!* 2 passata bis
!* si incrementa il contatore di 1
    do i = 1, Grid%nmax+1
       Icont_bp(i) = Icont_bp(i) + 1
    end do

!* 3 passata finale
!* si mettono le particelle al loro posto nel vettore NPartOrd
    do npi = 1,n_body_part
       NPartOrd_bp(Icont_bp(bp_arr(npi)%cell)+numpartincelgiaaposto(bp_arr(npi)%cell)) = npi
       numpartincelgiaaposto(bp_arr(npi)%cell) = numpartincelgiaaposto(bp_arr(npi)%cell) + 1
    end do
! 
 endif
!AA501b end

return
end subroutine OrdGrid1
!---split

