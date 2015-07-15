!cfile GetVarPart.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : GetVarPart
!
! Last updating : May 08, 2012
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to act on pl(0) where are pre-fixed x,y,z and define the
!                  particle variables
!
! Calling routine: CalcVarp
!
! Called routines: 
!
!************************************************************************************
!
subroutine GetVarPart (pglocal)
!* agisce sulla pl(0) di cui sono prefissati x,y,h
!* e definisce le variabili della particella
!
!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
!
!AA406 sub
integer(4)        :: nceli,ncel,npar,igridi,kgridi,jgridi,irang,krang,jrang,fw
!
integer(4)        :: mm,npj,irestocell
double precision  :: rijlocal,uni, pesoj,plocal,ro   
!
!.. Local Arrays ..
double precision, dimension (3) :: raglocal,vel
type (TyCtlPoint) :: pglocal

integer(4),      external :: CellIndices, CellNumber
double precision,external :: w
!
!.. Executable Statements ..
!
 nceli = pglocal%cella
 if (nceli == 0) return
 irestocell = CellIndices(nceli,igridi,jgridi,kgridi)
!
 uni    = zero
 plocal = zero
 ro     = zero
 vel(:) = zero
!
 npar = 0
 do jrang = jgridi-1,jgridi+1    ! ---- a  loop sulle 9 celle 
   do irang = igridi-1,igridi+1    ! ---- b  loop sulle 9 celle  
     do krang = kgridi-1,kgridi+1    ! ---- c  loop sulle 9 celle  
!
       ncel = CellNumber (irang,jrang,krang)
       if ( ncel == 0 ) cycle    ! cella fuori campo
!
       if (Icont(ncel+1) <= Icont(ncel)) cycle
       do mm = Icont(ncel),Icont(ncel+1)-1     ! loop sulle part di una cella  !20051230
         npj   = NPartOrd(mm)
!
         if ( pg(npj)%vel_type == "fix" ) cycle
!
         raglocal(:) = abs ( pglocal%coord(:)-pg(npj)%coord(:) )
         if ( raglocal(1) >= doubleh ) cycle
         if ( raglocal(2) >= doubleh ) cycle
         if ( raglocal(3) >= doubleh ) cycle
         rijlocal = raglocal(1)*raglocal(1) + raglocal(2)*raglocal(2) + raglocal(3)*raglocal(3)
         if ( rijlocal > doublesquareh ) cycle
         npar = npar+1
!
!AA406test
         rijlocal = dsqrt(rijlocal)       
!
         pesoj = pg(npj)%mass * w(rijlocal,Domain%h,Domain%coefke) / pg(npj)%dens
!
!.. calcolo den
         uni   = uni + pesoj
!
!.. calcolo vars
         plocal  = plocal  + pg(npj)%pres * pesoj        ! p locale
         ro = ro + pg(npj)%dens * pesoj                  ! ro VA FATTO SOLO SULLO STESSO MEZZO
         vel(:) = vel(:) + pg(npj)%vel(:) * pesoj        ! vx
!
       end do
!    
!AA406 start
       if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
! Loop over the neighbouring wall particles in the cell
          do fw = Icont_w(ncel),Icont_w(ncel+1)-1
             npj = NPartOrd_w(fw)
! Relative positions and distances
             raglocal(1:3) = pglocal%coord(:) - pg_w(npj)%coord(1:3)
             rijlocal = raglocal(1)*raglocal(1) + raglocal(2)*raglocal(2) + raglocal(3)*raglocal(3)
! Distance check
             if (rijlocal > doublesquareh) cycle
             rijlocal = dsqrt(rijlocal)
             pesoj = pg_w(npj)%mass * w(rijlocal,Domain%h,Domain%coefke) / pg_w(npj)%dens
             uni   = uni + pesoj
             plocal  = plocal  + pg_w(npj)%pres * pesoj
             ro = ro + pg_w(npj)%dens * pesoj                  
             vel(:) = vel(:) + pg_w(npj)%vel(:) * pesoj        
         end do 
      endif
!AA406 end
       
!
     end do
   end do
 end do
!
 pglocal%uni = uni
!!!! if ( npar > 0 ) then
! prova correzione per escludere punti che stanno al di fuori del campo
!
!AA501 sub
 if ( npar > 0 .and. uni >= 0.1) then
!
    pglocal%pres   = plocal/uni
    pglocal%dens   = ro/uni
    pglocal%vel(:) = vel(:)/uni
 end if
!
return
end subroutine GetVarPart
!---split

