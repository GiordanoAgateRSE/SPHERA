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
! Program unit: SubCalcPreIdro      
! Description: Hydrostatic pressure profiles (in case they are imposed as initial conditions).                 
!----------------------------------------------------------------------------------------------------------------------------------

subroutine SubCalcPreIdro
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
logical          :: foundcell
integer(4) :: idpel,Nz,pl_id,npi,iappo,i,j,k,numcell,ncelcorr,m,nnlocal,igridi
integer(4) :: jgridi,kgridi,nnsave
double precision :: ZQuotaMediumCorr,ZQuotaColonna,ZQuotaSecondMedium,affond1
double precision :: affond2,pl_quote,gravmod,coshor,senhor
character(len=lencard)  :: nomsub = "SubCalcPreIdro"
integer(4),external :: CellNumber,ParticleCellNumber,CellIndices
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
! To evaluate the absolute value of gravity and the unity vector aligned with 
! the vertical direction.
gravmod = Dsqrt(Domain%grav(1) * Domain%grav(1) + Domain%grav(2) *             &
          Domain%grav(2) + Domain%grav(3) * Domain%grav(3))
if (gravmod>zero) then
   coshor = - Domain%grav(3) / gravmod
   senhor = Domain%grav(1) / gravmod
   else
      coshor = zero
      senhor = zero
end if
pl_quote = max_negative_number
pl_id = 0
! To search for free level conditions, if present ("pl" type).
do npi=1,nag
   Nz = pg(npi)%izona
   if (partz(Nz)%pressure/="pl" ) cycle 
   if (pl_id==0) then
      pl_id = int(partz(Nz)%valp)
      else if (pl_id/=int(partz(Nz)%valp)) then
         call diagnostic (arg1=10,arg2=8,arg3=nomsub)
   end if
   pl_quote = max(pl_quote,pg(npi)%coord(3))
end do 
particle_loop: do npi=1,nag
! Initial conditions for pressure (from input file)
   Nz = pg(npi)%izona
   if (partz(Nz)%pressure=="pa" ) then  
      pg(npi)%pres = partz(Nz)%valp
      pg(npi)%dens = med(pg(npi)%imed)%den0
      cycle particle_loop
   end if   
! To detect the cell number of the current particle
   ncelcorr = ParticleCellNumber(pg(npi)%coord)
   iappo = CellIndices(ncelcorr,igridi,jgridi,kgridi)
   i = igridi
   j = jgridi
! To check if there is a cell including a medium interface in the column.
! If this is so, it evaluates the medium interface level.
   idpel = 0
   foundcell = .FALSE.
   ZQuotaMediumCorr = pg(npi)%coord(3)
   ZQuotaColonna = pg(npi)%coord(3)
   ZQuotaSecondMedium = pg(npi)%coord(3)
! Loop over the column cells above
   do k = kgridi,Grid%ncd(3)
      numcell = CellNumber(i,j,k)
      if (numcell==0) cycle
! Loop over the number of particles inside the cell
      do m=Icont(numcell),Icont(numcell+1)-1
         nnlocal = npartord(m)
! A particle of different medium is found, so the interface cell identifier 
! is set
         if ((med(pg(npi)%imed)%index/=med(pg(nnlocal)%imed)%index).and.       &
         (index(Med(pg(nnlocal)%imed)%tipo,"gas")==0)) then
! The interface cell is set
            if (.not.foundcell) then
               idpel = numcell
               foundcell = .TRUE.
               nnsave = nnlocal
            end if
! The minimum level in the interface cell for the medium different from the 
! current one is set
            if (numcell==idpel) then
               ZQuotaSecondMedium = min(ZQuotaSecondMedium,pg(nnlocal)%coord(3))
            end if
! To increase in any case the reference level of the column 
            ZQuotaColonna = max(ZQuotaColonna,pg(nnlocal)%coord(3))
            else
! To evaluate the reference levels: "ZQuotaMediumCorr" is the maximum level 
! of the same medium of the current particle, while "ZQuotaColonna" is the 
! maximum level of the other medium located above (maximum of the column of 
! cells)
               ZQuotaMediumCorr = max(ZQuotaMediumCorr,pg(nnlocal)%coord(3))
               ZQuotaColonna    = max(ZQuotaColonna,pg(nnlocal)%coord(3))
         end if
      end do
   end do
! To check if the current particles is inside the intermediate cell, but it is 
! of the upper medium type
   if (abs(ZQuotaMediumCorr-ZQuotaColonna)<xyz_tolerance) foundcell = .false.
! To set the reference pressure quote depending on the condition type
   if (partz(Nz)%pressure=="qp") ZQuotaColonna = partz(Nz)%valp
   if (partz(Nz)%pressure=="pl") ZQuotaColonna = pl_quote
! An upper medium interface cell has been found 
   if (foundcell) then
! The average quote of the medium interface in the cell is calculated
      ZQuotaMediumCorr  = (ZquotaMediumCorr + ZQuotaSecondMedium) * half
! particle pressure and density are evaluated, accounting
! for the medium interface
      affond1 = (ZQuotaColonna - ZQuotaMediumCorr) * coshor +                  &
                pg(nnsave)%coord(1) * senhor
      affond2 = (ZQuotaMediumCorr - pg(npi)%coord(3)) * coshor +               &
                pg(npi)%coord(1) * senhor 
      if (Domain%tipo == "bsph") then  
! To check this line 
         pg(npi)%pres = 0.d0 * (affond1 * med(pg(npi)%imed)%den0 * gravmod) *  &
                       (1.d0 - pg(npi)%coord(1) / 0.5925d0)
         else
            pg(npi)%pres = (affond1 * Med(pg(nnsave)%imed)%den0 + affond2 *    &
                            med(pg(npi)%imed)%den0) * gravmod
      endif
      else
         affond1 = (ZQuotaColonna - pg(npi)%coord(3)) * coshor +               &
                   pg(npi)%coord(1) * senhor 
         pg(npi)%pres = affond1 * med(pg(npi)%imed)%den0 * gravmod
         if (Domain%tipo == "bsph") then
! To check this line
            pg(npi)%pres = 0.d0 * (affond1 * med(pg(npi)%imed)%den0 * gravmod) &
                           * (1.d0 - pg(npi)%coord(1) / 0.5925d0)
            else
               pg(npi)%pres = (affond1 * med(pg(npi)%imed)%den0 * gravmod)
         end if
   endif
! Density
   pg(npi)%dens = (one + pg(npi)%pres/med(pg(npi)%imed)%eps) *                 &
                  med(pg(npi)%imed)%den0
   pg(npi)%dden = zero
! Diffusion model
   if (diffusione) then
      if ( pg(npi)%VolFra == VFmx) then
         pg(npi)%pres = ((ZQuotaColonna - ZQuotaMediumCorr) * (med(2)%den0 *   &
                        VFmn + med(1)%den0 * (one - VFmn)) + (ZQuotaMediumCorr &
                        - pg(npi)%coord(3)) * (med(2)%den0 * VFmx + med(1)%den0&
                        * (one - VFmx))) * gravmod 
         pg(npi)%rhoc = (one + pg(npi)%pres / med(pg(npi)%imed)%eps) *         &
                        med(pg(npi)%imed)%den0 
         pg(npi)%rhow = (one + pg(npi)%pres / med(1)%eps) * med(1)%den0
         pg(npi)%dens = VFmx * pg(npi)%rhoc + (one - VFmx) * pg(npi)%rhow
         else if (pg(npi)%VolFra==VFmn) then
            pg(npi)%pres = ((ZQuotaColonna - pg(npi)%coord(3)) * (med(2)%den0 *&
                           VFmn + med(1)%den0 * (one - VFmn))) * gravmod 
            pg(npi)%rhoc = (one + pg(npi)%pres / med(2)%eps) * med(2)%den0 
            pg(npi)%rhow = (one + pg(npi)%pres / med(pg(npi)%imed)%eps) *      &
                           med(pg(npi)%imed)%den0
            pg(npi)%dens = VFmn * pg(npi)%rhoc + (one - VFmn) * pg(npi)%rhow
      end if
   end if
enddo particle_loop
!------------------------
! Deallocations
!------------------------
return
end subroutine subCalcPreIdro

