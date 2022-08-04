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
! Program unit: SubCalcPreIdro      
! Description: To initialize pressure and density fields. Hydrostatic pressure 
!              profiles or unifrom pressure as initial conditions. Hydrostatic 
!              profiles are formally correct only in the presence of maximum 2 
!              fluid media along the vertical.            
!-------------------------------------------------------------------------------
subroutine SubCalcPreIdro
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: foundcell
integer(4) :: idpel,Nz,pl_id,npi,iappo,i,j,k,numcell,ncelcorr,m,nnlocal,igridi
integer(4) :: jgridi,kgridi,nnsave
double precision :: ZQuotaMediumCorr,ZQuotaColonna,ZQuotaSecondMedium,affond1
double precision :: affond2,pl_quote,gravmod,coshor,senhor
character(len=lencard)  :: nomsub = "SubCalcPreIdro"
integer(4),external :: CellNumber,ParticleCellNumber,CellIndices
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine EoS_barotropic_linear(k_bulk,rho_ref,p_ref,rho_in,p_in,rho_out,  &
      p_out)
      implicit none
      double precision,intent(in) :: k_bulk,rho_ref,p_ref
      double precision,intent(in),optional :: rho_in,p_in
      double precision,intent(out),optional :: rho_out,p_out
   end subroutine EoS_barotropic_linear
end interface
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
gravmod = dsqrt(Domain%grav(1) * Domain%grav(1) + Domain%grav(2) *             &
          Domain%grav(2) + Domain%grav(3) * Domain%grav(3))
if (gravmod>zero) then
   coshor = - Domain%grav(3) / gravmod
   senhor = Domain%grav(1) / gravmod
   else
      coshor = zero
      senhor = zero
endif
pl_quote = max_negative_number
pl_id = 0
! To search for free level conditions, if present ("pl" type).
do npi=1,nag
   Nz = pg(npi)%izona
   if (partz(Nz)%pressure/="pl") cycle
   if (pl_id==0) then
      pl_id = int(partz(Nz)%valp)
      elseif (pl_id/=int(partz(Nz)%valp)) then
         call diagnostic(arg1=10,arg2=8,arg3=nomsub)
   endif
   pl_quote = max(pl_quote,pg(npi)%coord(3))
enddo 
particle_loop: do npi=1,nag
! Initial conditions for pressure (from input file)
   Nz = pg(npi)%izona
   if (partz(Nz)%pressure=="pa") then  
      pg(npi)%pres = partz(Nz)%valp + Domain%prif
      pg(npi)%dens = med(pg(npi)%imed)%den0
      cycle particle_loop
   endif   
! To detect the cell number of the current particle
   ncelcorr = ParticleCellNumber(pg(npi)%coord)
   iappo = CellIndices(ncelcorr,igridi,jgridi,kgridi)
   i = igridi
   j = jgridi
! To check if there is a cell including a medium interface in the column.
! If this is so, it evaluates the medium interface level.
   idpel = 0
   foundcell = .false.
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
         if ((med(pg(npi)%imed)%index/=med(pg(nnlocal)%imed)%index)) then
! The interface cell is set
            if (.not.foundcell) then
               idpel = numcell
               foundcell = .true.
               nnsave = nnlocal
            endif
! The minimum level in the interface cell for the medium different from the 
! current one is set
            if (numcell==idpel) then
               ZQuotaSecondMedium = min(ZQuotaSecondMedium,pg(nnlocal)%coord(3))
            endif
! To increase in any case the reference level of the column 
            ZQuotaColonna = max(ZQuotaColonna,pg(nnlocal)%coord(3))
            else
! To evaluate the reference levels: "ZQuotaMediumCorr" is the maximum level 
! of the same medium of the current particle, while "ZQuotaColonna" is the 
! maximum level of the other medium located above (maximum of the column of 
! cells)
               ZQuotaMediumCorr = max(ZQuotaMediumCorr,pg(nnlocal)%coord(3))
               ZQuotaColonna    = max(ZQuotaColonna,pg(nnlocal)%coord(3))
         endif
      enddo
   enddo
! To check if the current particles is inside the intermediate cell, but it is 
! of the upper medium type
   if (abs(ZQuotaMediumCorr-ZQuotaColonna)<xyz_tolerance) foundcell = .false.
! To set the reference pressure quote depending on the condition type
   if (partz(Nz)%pressure=="qp") ZQuotaColonna = partz(Nz)%valp
   if (partz(Nz)%pressure=="pl") ZQuotaColonna = pl_quote + Domain%dx / 2.d0
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
      if (Domain%tipo=="bsph") then  
! To check this line
         pg(npi)%pres = 0.d0 * (affond1 * med(pg(npi)%imed)%den0 * gravmod) *  &
                       (1.d0 - pg(npi)%coord(1) / 0.5925d0)
         else
            pg(npi)%pres = (affond1 * Med(pg(nnsave)%imed)%den0 + affond2 *    &
                            med(pg(npi)%imed)%den0) * gravmod + Domain%prif
      endif
      else
         affond1 = (ZQuotaColonna - pg(npi)%coord(3)) * coshor +               &
                   pg(npi)%coord(1) * senhor 
         pg(npi)%pres = affond1 * med(pg(npi)%imed)%den0 * gravmod + Domain%prif
         if (Domain%tipo == "bsph") then
! To check this line
            pg(npi)%pres = 0.d0 * (affond1 * med(pg(npi)%imed)%den0 * gravmod) &
                           * (1.d0 - pg(npi)%coord(1) / 0.5925d0)
            else
               pg(npi)%pres = (affond1 * med(pg(npi)%imed)%den0 * gravmod) +   &
                              Domain%prif
         endif
   endif
! EoS inverse
   call EoS_barotropic_linear(Med(pg(npi)%imed)%eps,Med(pg(npi)%imed)%den0,    &
      Domain%prif,p_in=pg(npi)%pres,rho_out=pg(npi)%dens)
   pg(npi)%dden = zero
enddo particle_loop
!------------------------
! Deallocations
!------------------------
return
end subroutine subCalcPreIdro
