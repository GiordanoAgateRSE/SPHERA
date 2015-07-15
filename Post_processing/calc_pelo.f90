!cfile calc_pelo.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : calc_pelo
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to write results for free surface
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
subroutine calc_pelo
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4)       :: i,j, ncel,npartcel,mm,jj  !!!!!,icord,n 
!!!!!double precision :: ro1, ro2, romezzi
!!!!!logical          :: trovato
!
!.. Local Arrays ..
integer(4),      dimension(1)  :: minpos1,minpos2
double precision,dimension(3)  :: ragtemp
double precision,dimension(3,nlines)  :: pelolib
integer(4),      dimension(:),allocatable :: PartCelnum
double precision,dimension(:),allocatable :: PartCel
 
!
!.. Executable Statements ..
!
!!!!! pelolib = zero
!!!!! romezzi = Med(Domain%ipllb_md)%den0 * half
!!!!! trovato = .false.
!!!!!
!!!!! do i = 1,nlines
!!!!!
!!!!!    POINTS_LOOP: do j = control_lines(i)%Icont(1)+1,control_lines(i)%Icont(2)
!!!!!!       do n = 1,ncord
!!!!!!          icord = icoordp(n,ncord-1)
!!!!!!          pelolib(n,i) = control_points(j-1)%coord(icord)
!!!!!!       end do
!!!!!!       ro1 = control_points(j-1)%dens * control_points(j-1)%uni
!!!!!!       ro2 = control_points(j  )%dens * control_points(j  )%uni
!!!!!       ro1 = control_points(j-1)%dens * half
!!!!!       ro2 = control_points(j  )%dens * half
!!!!!
!!!!!       if ( ro2 < romezzi ) then
!!!!!          trovato = .true.
!!!!!          do n = 1,ncord
!!!!!             icord = icoordp(n,ncord-1)
!!!!!             pelolib(n,i) = control_points(j-1)%coord(icord) + &
!!!!!                            ( romezzi - ro1 ) / ( ro2 - ro1 ) * &
!!!!!                            ( control_points(j)%coord(icord) - control_points(j-1)%coord(icord) )
!!!!!          end do
!!!!!          exit POINTS_LOOP
!!!!!       end if
!!!!!    end do POINTS_LOOP
!!!!! end do
!!!!!
!!!!! if (trovato) write (nplb,'(30g14.7)') tempo,pelolib

 allocate (PartCelnum(NMAXPARTJ), PartCel(NMAXPARTJ))
!.. loop su tutte le linee
 Pelolib = 0
 do i = 1,nlines
!.. loop sui punti della linea
   POINTS_LOOP: do j = control_lines(i)%Icont(1)+1,control_lines(i)%Icont(2)
     ncel = control_points(j)%cella
     if (ncel == 0) cycle    ! cella fuori campo
!.. se la cella a cui appartiene il punto non contiene particelle esco
     if (Icont(ncel+1) <= Icont(ncel)) cycle
!.. altrimenti trovo le due particelle piu' vicine al punto e faccio la media della posizione
     nPartCel = 0
     PartCelnum = 0
     PartCel = 99999
     do mm = Icont(ncel),Icont(ncel+1)-1  ! +++ loop sulle part di una cella
       jj = NPartOrd(mm)
!.. calcola le conponenti deltaX, deltaY e deltaZ della distanza tra il punto e la particella
       ragtemp(1:3) = control_points(j)%coord(1:3) - pg(jj)%coord(1:3)
!.. calcola la distanza effettiva (considera direttamente il quadrato per aumentare l'accuratezza)
       npartcel = npartcel + 1
       PartCelnum(npartcel) = jj
       PartCel(npartcel) = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
     end do  ! +++ loop sulle part di una cella
!.. cerco le particelle piu' vicine al punto
     minpos1 = minloc(PartCel,1)
     PartCel(minpos1) = 9999.
     minpos2 = minloc(PartCel,1)
 !.. posizione del pelo libero come media della posizione delle due particelle piu' vicine al punto
     pelolib(1:3,i) = (pg(PartCelnum(minpos1(1)))%coord(1:3) + pg(PartCelnum(minpos2(1)))%coord(1:3)) * half
   end do POINTS_LOOP
 end do
 write (nplb,'(30g14.7)') tempo,pelolib
 deallocate (PartCelnum,PartCel)

return
end subroutine calc_pelo
!---split

