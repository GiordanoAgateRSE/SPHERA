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
! Program unit: MohrC
! Description: Mohr-Coulomb 2D erosion criterion (Manenti et al., 2012, JHE). 
!              Shield erosion criterion works better (Manenti et al., 2012, JHE).                   
!----------------------------------------------------------------------------------------------------------------------------------

subroutine MohrC 
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
integer(4) :: npi,j,k,m,intpl_id,imed,pl_imed,intliq_id,intsol_id,ncelcorr     
integer(4) :: igridi,jgridi,kgridi,iappo,numcell,nnlocal,flag,nsp,npartint,npj
double precision :: PartLiq_pres,peloloc,interf_liq,interf_sol,preidro,pretot
double precision :: preeff,Velocity2,appo1,secinv,coeff1,coeff2,mu,mumax
character(len=lencard)  :: nomsub = "Mohr-Coulomb"
integer(4),external :: ParticleCellNumber, CellIndices, CellNumber
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
!$omp parallel do default(none)                                                &
!$omp private(npi,imed,ncelcorr,iappo,igridi,jgridi,kgridi,flag,nsp,appo1,j)   &
!$omp private(npartint,npj,k,numcell,nnlocal,m,intpl_id,peloloc,pl_imed)       &
!$omp private(interf_liq,intliq_id,interf_sol,intsol_id,secinv,PartLiq_pres)   &
!$omp private(preidro,pretot,preeff,coeff1,coeff2,mu,mumax,Velocity2)          &
!$omp shared(ind_interfaces,nomsub,nout,nscr,nag,pg,Med,Domain,Grid,Icont)     &
!$omp shared(npartord,nPartIntorno,PartIntorno,NMAXPARTJ,diffusione,esplosione)&
!$omp shared(it_corrente)
do npi = 1,nag
   if ((pg(npi)%cella==0).or.(pg(npi)%vel_type/="std")) cycle
! Check motion status and viscosity computation 
   imed = pg(npi)%imed
   if (index(Med(imed)%tipo,"granular")<=0) cycle
! In case of explosion, erosion criterion is not active. 
   if (esplosione) then
      Velocity2 = pg(npi)%vel(1) * pg(npi)%vel(1) + pg(npi)%vel(2) *           &
                  pg(npi)%vel(2) + pg(npi)%vel(3) * pg(npi)%vel(3)
      if (Velocity2>1.0e-3) cycle
   end if
   ncelcorr = ParticleCellNumber(pg(npi)%coord)
   iappo = CellIndices(ncelcorr,igridi,jgridi,kgridi)
   intpl_id = ind_interfaces(igridi,jgridi,1)
   intliq_id = ind_interfaces(igridi,jgridi,2)
   intsol_id = ind_interfaces(igridi,jgridi,3)
   if (intliq_id==0) then
! Detection free surface  
      interf_liq = max_positive_number
      intliq_id = 0
! Loop over the on-going cell and the cells above 
      do k=1,Grid%ncd(3)
         numcell = CellNumber(igridi,jgridi,k)
         if (numcell==0) cycle
! Loop over the particles of the cell "numcell"
         do m=Icont(numcell),Icont(numcell+1)-1
            nnlocal = npartord(m)
            if (index(Med(pg(nnlocal)%imed)%tipo,"liquid")>0) then
! Detection liquid/granular interface  
               if (pg(nnlocal)%Coord(3)<interf_liq) then
                  interf_liq = pg(nnlocal)%Coord(3)
                  intliq_id = nnlocal
               end if
            end if
         end do
      end do
   end if
   if (intsol_id==0) then
! Detection solid interface 
      interf_sol = max_negative_number
      intsol_id = 0
! Loop over the on-going cell and the cells above 
      do k=1,Grid%ncd(3)
         numcell = CellNumber(igridi,jgridi,k)
         if (numcell==0) cycle
! Loop over the particles of the cell "numcell"
         do m = Icont(numcell),Icont(numcell+1)-1
            nnlocal = npartord(m)
            if (index(Med(pg(nnlocal)%imed)%tipo,"granular")>0) then
! Detection liquid/granular interface (granular particle)
               if (pg(nnlocal)%Coord(3)>interf_sol) then
                  interf_sol = pg(nnlocal)%Coord(3)
                  intsol_id = nnlocal
               endif
            endif
         enddo
      end do
   end if
   if ((intliq_id==0).and.(intsol_id==0)) then
      call diagnostic (arg1=11,arg2=2,arg3=nomsub)
   end if
   if (intpl_id==0) then
! Detection of free surface 
      peloloc = max_negative_number
      intpl_id = 0
! Loop over the on-going cell and the cells above 
      do k = 1,Grid%ncd(3)
         numcell = CellNumber(igridi,jgridi,k)
         if (numcell==0) cycle
! Loop over the particles of the cell "numcell"
         do m = Icont(numcell),Icont(numcell+1)-1
            nnlocal = npartord(m)
            if (index(Med(pg(nnlocal)%imed)%tipo,"liquid")>0) then
! Detection liquid/granular interface (liquid particle)
               if (pg(nnlocal)%Coord(3)>peloloc) then
                  peloloc = pg(nnlocal)%Coord(3)
                  intpl_id = nnlocal
               endif
            endif
         enddo
      enddo
   end if
   if (intliq_id==0) intliq_id  = intsol_id
   interf_liq = pg(intliq_id)%coord(3)
   if (intsol_id==0) intsol_id = intliq_id
   interf_sol = pg(intsol_id)%coord(3)
   if (intpl_id==0) intpl_id = intliq_id
   peloloc = pg(intpl_id)%coord(3)
   pl_imed = pg(intpl_id)%imed
! Apparent viscosity
   secinv = two * pg(npi)%secinv
! secinv depending on the viscosity at the free surface 
   secinv = secinv * Med(imed)%mumx / (Med(pl_imed)%visc*med(pl_imed)%den0)
   preidro = - Domain%grav(3) * med(pl_imed)%den0 * (peloloc -                 &
             pg(npi)%Coord(3))
   Velocity2 = pg(intliq_id)%vel(1) * pg(intliq_id)%vel(1) +                   & 
             pg(intliq_id)%vel(2) * pg(intliq_id)%vel(2) +                     &
             pg(intliq_id)%vel(3)*pg(intliq_id)%vel(3)
   PartLiq_pres = (interf_liq - interf_sol) * med(pg(intliq_id)%imed)%den0 *   &
                  ( - Domain%grav(3)) + pg(intliq_id)%pres + Velocity2 * half  &
                  * med(pg(intliq_id)%imed)%den0
   pretot  = PartLiq_pres + med(pg(intsol_id)%imed)%den0 *                     &
             ( - Domain%grav(3)) * (interf_sol - pg(npi)%Coord(3))
   preeff  = pretot - preidro
   coeff1  = cos (Med(imed)%phi)
   coeff2  = sin (Med(imed)%phi)
   mu = (Med(imed)%coes * coeff1 + preeff * coeff2) / (secinv + 0.00001d0)
   mumax = Med(imed)%mumx
! "flag" is 1 if a 'liquid' particle belongs to the kernel support of the 
! on-going granular particle 
   flag = 0  
   if (mu<mumax) then 
      Nsp = nPartIntorno(npi) 
! Searching "liquid" neighbours to detect "solid" particles close to the 
! fluid-solid interface 
      appo1 = peloloc - interf_sol
      do j=1,Nsp
         npartint = (npi - 1)* NMAXPARTJ + j
         npj = PartIntorno(npartint)
! To include granular particles with motion status equal to "flu"
         if ((index(Med(pg(npj)%imed)%tipo,"liquid")>0).or.                    &
            ((index(Med(pg(npj)%imed)%tipo,"granular")>0).and.                 &
            (pg(npi)%state=="flu"))) then
            if (appo1>=Domain%dd) then
               flag = 1
               exit
            end if
         end if
      enddo
! Search completed                          
      if ((flag==0).and.(it_corrente>Med(imed)%NIterSol)) then
         if (pg(npi)%CloseBcOut==0) then
            pg(npi)%state = "sol"
            pg(npi)%vel = zero
            pg(npi)%var = zero
! Density consistent with hydrostatic pressure 
            if (.not.diffusione) pg(npi)%dens = med(imed)%den0 + (pretot /     &
                                                (Med(imed)%celerita *          &
                                                Med(imed)%celerita))
         end if
         else if (flag==1) then
            pg(npi)%state = "flu"
      end if
      else
         if (pg(npi)%CloseBcOut==0) then
            pg(npi)%state = "sol"
            pg(npi)%vel = zero
            pg(npi)%var = zero
! Density consistent with hydrostatic pressure 
            if (.not. diffusione) pg(npi)%dens = med(imed)%den0 + (pretot /    &
                                                 (Med(imed)%celerita *         &
                                                 Med(imed)%celerita))
         end if
   end if
end do
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
return
end subroutine MohrC

