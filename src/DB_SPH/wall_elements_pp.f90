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
! Program unit: wall_elements_pp
! Description: Smoothing wall element values for post-processing. 
!              Post-processing the wall surface element values (provided a
!              selected region). Post-processing the hydrodynamic normal
!              force on DBSPH surface elements (provided a selected region).
!              Post-processing the wall surface element values (provided
!              selected element IDs).               
!-------------------------------------------------------------------------------
subroutine wall_elements_pp
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: igridi,kgridi,jgridi,irang,krang,jrang,ncelj,jgrid1,jgrid2,i
integer(4) :: npi,npj,npartint,index_rij_su_h,irestocell,ww
double precision :: rij_su_h,ke_coef,rij_su_h_quad,rijtemp,rijtemp2,wu,denom
double precision :: gradmod,Fx
double precision :: ragtemp(3)
double precision, dimension(:), allocatable :: pres_wpp,den_wpp
character(255)   :: nomefilectl_wall
integer(4),external :: CellIndices, CellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
allocate(pres_wpp(DBSPH%n_w))
allocate(den_wpp(DBSPH%n_w))
!------------------------
! Initializations
!------------------------
Fx = 0.d0
ke_coef = Domain%coefke / Domain%h
pres_wpp = zero
den_wpp = zero
!------------------------
! Statements
!------------------------
! Loop over the wall elements
!$omp parallel do default(none)                                                &
!$omp shared(DBSPH,pg_w,ncord,Icont_w,NPartOrd_w,pres_wpp,Domain,den_wpp)      &
!$omp shared(doublesquareh,squareh)                                            &
!$omp private(npi,irestocell,igridi,jgridi,jgrid1,jgrid2,kgridi,jrang,irang)   &
!$omp private(krang,ncelj,ww,npj,rijtemp,ragtemp,rijtemp2,rij_su_h)            &
!$omp private(rij_su_h_quad,index_rij_su_h,wu)
do npi = 1,DBSPH%n_w
   if (pg_w(npi)%cella==0) cycle
   irestocell = CellIndices (pg_w(npi)%cella,igridi,jgridi,kgridi)
! Indices for grid searching in y-direction 
   jgrid1 = jgridi - (ncord - 2)
   jgrid2 = jgridi + (ncord - 2)
! Loop over 9 cells (search along x-direction)
   loop_jrang: do jrang = jgrid1,jgrid2     
! Loop over 1/9 cells (search along y-direction)
      loop_irang: do irang = igridi-1,igridi+1      
! Loop over 9 cells (search along z-direction)  
        loop_krang: do krang = kgridi-1,kgridi+1    
           ncelj = CellNumber(irang,jrang,krang)
! Cell out of domain 
           if (ncelj==0) cycle    
           if (Icont_w(ncelj+1)<=Icont_w(ncelj)) cycle
! Loop over the neighbouring wall particles in the cell
           do ww = Icont_w(ncelj),Icont_w(ncelj+1)-1
              npj = NPartOrd_w(ww)
! Relative positions and distances
              ragtemp(1:3) = pg_w(npi)%coord(1:3) - pg_w(npj)%coord(1:3)
              rijtemp = ragtemp(1) * ragtemp(1) + ragtemp(2) * ragtemp(2) +    &
                        ragtemp(3) * ragtemp(3)
! Distance check
              if (rijtemp > doublesquareh) cycle
              rijtemp2 = rijtemp
              rij_su_h = Dsqrt(rijtemp) / Domain%h
              rij_su_h_quad = rijtemp2 / squareh
              index_rij_su_h = int(rij_su_h)
! Kernel computation
              if (index_rij_su_h>=2) cycle
              wu = 0.666666667d0 + rij_su_h_quad * (rij_su_h * half - one)
              if (index_rij_su_h>0) wu = (two - rij_su_h) * (two - rij_su_h)   &
                                         * (two - rij_su_h) * 0.166666666667d0
              if (pg_w(npj)%wet==1) then
                 pres_wpp(npi) = pres_wpp(npi) + wu * Domain%coefke *          &
                                 pg_w(npj)%pres * pg_w(npj)%weight
                 den_wpp(npi)  = den_wpp(npi) +  wu * Domain%coefke *          &
                                 pg_w(npj)%weight
              endif
           end do
        end do loop_krang       
      end do loop_irang       
   end do loop_jrang   
end do 
!$omp end parallel do
! Loop over the wall particles to update the pressure values for post-processing
!$omp parallel do default(none) shared(DBSPH,den_wpp,pres_wpp) private(npi)
do npi=1,DBSPH%n_w
   if (den_wpp(npi)>0.d0) pres_wpp(npi) = pres_wpp(npi) / den_wpp(npi) 
end do
!$omp end parallel do
! Writing the pressure force (post-processing) (provided a region)
if (DBSPH%n_monitor_regions==1) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),      &
                                          '_wall_Fx_',on_going_time_step,".txt"
   open (unit_dbsph_Fx,file=nomefilectl_wall,status="unknown",form="formatted")
   write (unit_dbsph_Fx,*) "Force "
   write (unit_dbsph_Fx,'(1x,2(a,10x))') " Time(s)","Fx(kgm/s^2)"
   call flush(unit_dbsph_Fx)
   write (unit_dbsph_Fx,*) "Force at boundaries (DBSPH)"
   do npi=1,DBSPH%n_w
      if ((pg_w(npi)%coord(1)>=DBSPH%monitor_region(1)).and.                   &
          (pg_w(npi)%coord(1)<=DBSPH%monitor_region(2)).and.                   &
          (pg_w(npi)%coord(2)>=DBSPH%monitor_region(3)).and.                   &
          (pg_w(npi)%coord(2)<=DBSPH%monitor_region(4)).and.                   &
          (pg_w(npi)%coord(3)>=DBSPH%monitor_region(5)).and.                   &
          (pg_w(npi)%coord(3)<=DBSPH%monitor_region(6))) then
         if (pg_w(npi)%normal(1)/=0.d0) Fx = Fx + pres_wpp(npi) *              &
                                           pg_w(npi)%normal(1)*pg_w(npi)%weight
      endif
   end do
   write (unit_dbsph_Fx,'(2(1x,g14.7))') simulation_time,Fx
   close (unit_dbsph_Fx)
endif
! Writing the wall element pressure values derived from post-processing 
! (provided a region)
if (DBSPH%n_monitor_regions==1) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),      &
      '_wall_region_',on_going_time_step,".txt"
   open (unit_dbsph_se_reg,file=nomefilectl_wall,status="unknown",             &
      form="formatted")
   write (unit_dbsph_se_reg,*) "Wall element values "
   write (unit_dbsph_se_reg,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)')            &
      " Time","Iter","ID","X Coord","Y Coord","Z Coord","X Velocity",          &
      "Y Velocity","Z Velocity"," Pressure"," Pressure_pp","wet"
   call flush(unit_dbsph_se_reg)
   write (unit_dbsph_se_reg,*) "wall_elements"
   do npi=1,DBSPH%n_w
      if ((pg_w(npi)%coord(1)>=DBSPH%monitor_region(1)).and.                   &
          (pg_w(npi)%coord(1)<=DBSPH%monitor_region(2)).and.                   &
          (pg_w(npi)%coord(2)>=DBSPH%monitor_region(3)).and.                   &
          (pg_w(npi)%coord(2)<=DBSPH%monitor_region(4)).and.                   &
          (pg_w(npi)%coord(3)>=DBSPH%monitor_region(5)).and.                   &
          (pg_w(npi)%coord(3)<=DBSPH%monitor_region(6))) then
         write (unit_dbsph_se_reg,'(g14.7,2(i14),8(1x,g14.7),i3)')             &
            simulation_time,on_going_time_step,i,pg_w(npi)%coord(:),           &
            pg_w(npi)%vel(:),pg_w(npi)%pres,pres_wpp(npi),pg_w(npi)%wet 
      endif
   end do
   close (unit_dbsph_se_reg)
endif
! Writing the wall element pressure values derived from post-processing 
! (provided the element IDs)
if (DBSPH%n_monitor_points>0) then
   write(nomefilectl_wall,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),      &
      '_wall_IDs_',on_going_time_step,".txt"
   open (unit_dbsph_se_ID,file=nomefilectl_wall,status="unknown",              &
      form="formatted")
   write (unit_dbsph_se_ID,*) "Wall element values "
   write (unit_dbsph_se_ID,'(1x,2(a,10x),3(a,8x),3(a,5x),a,7x,a)')             &
      " Time","Iter","ID","X Coord","Y Coord","Z Coord","X Velocity",          &
      "Y Velocity","Z Velocity"," Pressure"," Pressure_pp","wet"
   call flush(unit_dbsph_se_ID)
   write (unit_dbsph_se_ID,*) "wall_elements"
   do i=1,DBSPH%n_monitor_points
      write (unit_dbsph_se_ID,'(g14.7,2(i14),8(1x,g14.7),i3)') simulation_time,&
         on_going_time_step,i,pg_w(DBSPH%monitor_IDs(i))%coord(:),             &
         pg_w(DBSPH%monitor_IDs(i))%vel(:),pg_w(DBSPH%monitor_IDs(i))%pres,    &
         pres_wpp(DBSPH%monitor_IDs(i)),pg_w(DBSPH%monitor_IDs(i))%wet           
   end do   
   close (unit_dbsph_se_ID)
endif
!------------------------
! Deallocations
!------------------------
deallocate(pres_wpp)
deallocate(den_wpp)
return
end subroutine wall_elements_pp

