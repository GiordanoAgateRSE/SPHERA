!AA601 sub the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! File name     : wall_elements_pp
!
! Creation      : Amicarelli A.; 30Nov11
! Update        : Amicarelli A.; 26Jan15 
!                 Integration with no hard-coding; parallelization; 
!                 post-processing the wall surface element values (provided a selected region)
!                 post-processing the hydrodynamic normal force on DBSPH surface elements (provided a selected region)
!                 post-processing the wall surface element values (provided selected element IDs)
!
!************************************************************************************
! Module purpose : Smoothing wall element values for post-processing; initializing
!
! Calling routine: Memo_Ctl
!
!************************************************************************************

subroutine wall_elements_pp

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
integer(4)       :: igridi,kgridi,jgridi,irang,krang,jrang,ncelj,jgrid1,jgrid2,i
integer(4)       :: npi,npj,npartint,index_rij_su_h,irestocell,ww
double precision :: rij_su_h,ke_coef,rij_su_h_quad,rijtemp,rijtemp2,wu,denom,gradmod
double precision :: ragtemp(3)
double precision, dimension(:), allocatable :: pres_wpp,den_wpp
character(255)   :: nomefilectl_wall
integer(4),external :: CellIndices, CellNumber
double precision :: Fx

! Initializations
Fx = 0.
ke_coef = Domain%coefke / Domain%h
allocate (pres_wpp(DBSPH%n_w))
allocate (den_wpp(DBSPH%n_w))
pres_wpp = zero
den_wpp = zero

! Statements
! Loop over the wall elements
!$omp parallel do default(none) &
!$omp shared(DBSPH,pg_w,ncord,Icont_w,NPartOrd_w,pres_wpp,Domain,den_wpp,doublesquareh,squareh) &
!$omp private(npi,irestocell,igridi,jgridi,jgrid1,jgrid2,kgridi,jrang,irang,krang,ncelj,ww,npj,rijtemp,ragtemp,rijtemp2,rij_su_h,rij_su_h_quad,index_rij_su_h,wu)
do npi = 1,DBSPH%n_w
   if (pg_w(npi)%cella == 0) cycle
   irestocell = CellIndices (pg_w(npi)%cella,igridi,jgridi,kgridi)
! Indici inizio-fine scansione Y per casi 2d o 3d
   jgrid1 = jgridi - (ncord-2)
   jgrid2 = jgridi + (ncord-2)
   loop_jrang: do jrang = jgrid1,jgrid2    ! ---- a  loop sulle 9 celle 
      loop_irang: do irang = igridi-1,igridi+1    ! ---- b  loop sulle 9 celle  
        loop_krang: do krang = kgridi-1,kgridi+1    ! ---- c  loop sulle 9 celle  
           ncelj = CellNumber (irang,jrang,krang)
           if (ncelj == 0) cycle    ! cella fuori campo
           if (Icont_w(ncelj+1) <= Icont_w(ncelj)) cycle
! Loop over the neighbouring wall particles in the cell
           do ww = Icont_w(ncelj),Icont_w(ncelj+1)-1
              npj = NPartOrd_w(ww)
! Relative positions and distances
              ragtemp(1:3) = pg_w(npi)%coord(1:3) - pg_w(npj)%coord(1:3)
              rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) + ragtemp(3)*ragtemp(3)
! Distance check
              if (rijtemp > doublesquareh) cycle
              rijtemp2 = rijtemp
              rij_su_h = Dsqrt(rijtemp) / Domain%h
              rij_su_h_quad = rijtemp2 / squareh
              index_rij_su_h = int(rij_su_h)
! Kernel computation
              if (index_rij_su_h >= 2) cycle
              wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
              if (index_rij_su_h > 0) wu = (two - rij_su_h) * (two - rij_su_h) * (two - rij_su_h) * 0.166666666667d0
              if (pg_w(npj)%wet==1) then
                 pres_wpp(npi) = pres_wpp(npi) + wu * Domain%coefke * pg_w(npj)%pres * pg_w(npj)%weight
                 den_wpp(npi)  = den_wpp(npi) +  wu * Domain%coefke * pg_w(npj)%weight
              endif
           end do
        end do loop_krang   ! ---- c  loop sulle 9 celle    
      end do loop_irang   ! ---- b  loop sulle 9 celle    
   end do loop_jrang   ! ---- a  loop sulle 9 celle
end do 
!$omp end parallel do
! Loop over the wall particles to update the pressure values for post-processing
!$omp parallel do default(none) shared(DBSPH,den_wpp,pres_wpp) private(npi)
do npi=1,DBSPH%n_w
   if (den_wpp(npi) > 0.) pres_wpp(npi) = pres_wpp(npi) / den_wpp(npi) 
end do
!$omp end parallel do
! Writing the pressure force (post-processing) (provided a region)
if (DBSPH%n_monitor_regions==1) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_wall_Fx_',it_corrente,".txt"
   open (unit_dbsph_Fx, file=nomefilectl_wall, status="unknown", form="formatted" )
   write (unit_dbsph_Fx,*) "Force "
   write (unit_dbsph_Fx,'(1x,2(a,10x))') " Time(s)","Fx(kgm/s^2)"
   call flush(unit_dbsph_Fx)
   write (unit_dbsph_Fx,*) "Force at boundaries (DBSPH)"
   do npi=1,DBSPH%n_w
      if ( (pg_w(npi)%coord(1)>=DBSPH%monitor_region(1)) .and. (pg_w(npi)%coord(1)<=DBSPH%monitor_region(2)) .and. &
           (pg_w(npi)%coord(2)>=DBSPH%monitor_region(3)) .and. (pg_w(npi)%coord(2)<=DBSPH%monitor_region(4)) .and. &
           (pg_w(npi)%coord(3)>=DBSPH%monitor_region(5)) .and. (pg_w(npi)%coord(3)<=DBSPH%monitor_region(6))        ) then
         if (pg_w(npi)%normal(1)/=0.) Fx=Fx+pres_wpp(npi)*pg_w(npi)%normal(1)*pg_w(npi)%weight
      endif
   end do
   write (unit_dbsph_Fx,'(2(1x,g14.7))') tempo,Fx
   close (unit_dbsph_Fx)
endif
! Writing the wall element pressure values derived from post-processing (provided a region)
if (DBSPH%n_monitor_regions==1) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_wall_region_',it_corrente,".txt"
   open (unit_dbsph_se_reg, file=nomefilectl_wall, status="unknown", form="formatted" )
   write (unit_dbsph_se_reg,*) "Wall element values "
   write (unit_dbsph_se_reg,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
         " Time","Iter","ID","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure"," Pressure_pp","wet"
   call flush(unit_dbsph_se_reg)
   write (unit_dbsph_se_reg,*) "wall_elements"
   do npi=1,DBSPH%n_w
      if ( (pg_w(npi)%coord(1)>=DBSPH%monitor_region(1)) .and. (pg_w(npi)%coord(1)<=DBSPH%monitor_region(2)) .and. &
           (pg_w(npi)%coord(2)>=DBSPH%monitor_region(3)) .and. (pg_w(npi)%coord(2)<=DBSPH%monitor_region(4)) .and. &
           (pg_w(npi)%coord(3)>=DBSPH%monitor_region(5)) .and. (pg_w(npi)%coord(3)<=DBSPH%monitor_region(6))        ) then
         write (unit_dbsph_se_reg,'(g14.7,2(i14),8(1x,g14.7),i3)') tempo,it_corrente,i,pg_w(npi)%coord(:),pg_w(npi)%vel(:), &
                                                                   pg_w(npi)%pres,pres_wpp(npi),pg_w(npi)%wet 
      endif
   end do
   close (unit_dbsph_se_reg)
endif
! Writing the wall element pressure values derived from post-processing (provided the element IDs)
if (DBSPH%n_monitor_points>0) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_wall_IDs_',it_corrente,".txt"
   open (unit_dbsph_se_ID, file=nomefilectl_wall, status="unknown", form="formatted" )
   write (unit_dbsph_se_ID,*) "Wall element values "
   write (unit_dbsph_se_ID,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)') &
      " Time","Iter","ID","X Coord","Y Coord","Z Coord","X Velocity","Y Velocity","Z Velocity"," Pressure"," Pressure_pp","wet"
   call flush(unit_dbsph_se_ID)
   write (unit_dbsph_se_ID,*) "wall_elements"
   do i=1,DBSPH%n_monitor_points
      write (unit_dbsph_se_ID,'(g14.7,2(i14),8(1x,g14.7),i3)') tempo,it_corrente,i,pg_w(DBSPH%monitor_IDs(i))%coord(:),pg_w(DBSPH%monitor_IDs(i))%vel(:), &
                                                            pg_w(DBSPH%monitor_IDs(i))%pres,pres_wpp(DBSPH%monitor_IDs(i)),pg_w(DBSPH%monitor_IDs(i))%wet           
   end do   
   close (unit_dbsph_se_ID)
endif

! Deallocations
deallocate (pres_wpp)
deallocate (den_wpp)

return
end subroutine wall_elements_pp
!---split

