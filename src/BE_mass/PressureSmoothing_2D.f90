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
! Program unit: PressureSmoothing_2D
! Description: Partial smoothing for pressure (Di Monaco et al., 2011), also 
!              with DB-SPH boundary treatment scheme. This subroutine is not 
!              just the formal 2D extension of "PressureSmoothing_3D.f90": 
!              differences are appreciable.
!-------------------------------------------------------------------------------
#ifdef SPACE_2D
subroutine PressureSmoothing_2D
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: Ncbs,IntNcbs,npi,npj,mati,nsp,icbs,iside,sidestr,izonelocal     
integer(4) :: ibdt,ibdp,j,i,ii,npartint
double precision :: ro0i,p0i,pi,sompW,pesoj,TetaP1,VIntWdV_FT,VIntWdV_SO      
double precision :: VIntWdV_OSB,VIntWdV_OSP,press_so, press_osb, vin
double precision :: Appunity,smoothpi,IntWdV
integer(4),dimension(1:PLANEDIM) :: acix
double precision,dimension(1:2) :: IntDpWdV
double precision,dimension(1:PLANEDIM) :: IntLocXY,sss,nnn,massforce
#ifdef SOLID_BODIES
double precision,dimension(:),allocatable :: sompW_vec,AppUnity_vec
#endif
!="a"(SPHERA),"b"(FreeSurf),"c"
character(1),parameter :: SmoothVersion = "b"
character(4) :: strtype
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
#ifdef SOLID_BODIES
   allocate(sompW_vec(nag))
   allocate(AppUnity_vec(nag))
   sompW_vec = zero
   AppUnity_vec = zero
#endif
acix(1) = 1
acix(2) = 3
!------------------------
! Statements
!------------------------
#ifdef SOLID_BODIES
! Body particle contributions to pressure smoothing
   call start_and_stop(2,19)
   call body_to_smoothing_pres(sompW_vec,AppUnity_vec)
   call start_and_stop(3,19)
#endif
!$omp parallel do default(none)                                                &
!$omp shared(nag,Pg,Med,Tratto,Partz,Domain,nPartIntorno,PartIntorno,NMAXPARTJ)&
!$omp shared(PartKernel,kernel_fw,BoundarySide,BoundaryDataPointer)            &
!$omp shared(BoundaryDataTab,acix,dt,indarrayFlu,Array_Flu,input_any_t)        &
#ifdef SOLID_BODIES
!$omp shared(n_bodies,sompW_vec,AppUnity_vec)                                  &
#endif
!$omp private(npi,ii,Appunity,TetaP1,Ncbs,IntNcbs,ibdt,icbs,ibdp,iside)        &
!$omp private(sidestr,Nsp,mati,ro0i,p0i,pi,SompW,j,npartint,npj,pesoj,smoothpi)&
!$omp private(VIntWdV_FT,VIntWdV_SO,VIntWdV_OSB,VIntWdV_OSP,press_so,press_osb)&
!$omp private(IntLocXY,strtype,sss,nnn,massforce,IntWdV,IntDpWdV,izonelocal)   &
!$omp private(vin)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
! Excluding particles close to the face with conditions "flow", "velo" and "sour"
   if (pg(npi)%koddens==0) then 
      Nsp = nPartIntorno(npi)
      if (Nsp>0) then
         Appunity = zero
         mati = pg(npi)%imed
         ro0i = Med(mati)%den0
         p0i = Domain%prif
         pi = pg(npi)%pres
         sompW = zero
         do j=1,Nsp
            npartint = (npi - 1) * NMAXPARTJ + j
            npj = PartIntorno(npartint)
            pesoj = pg(npj)%mass * PartKernel(4,npartint) / pg(npj)%dens
            Appunity = Appunity + pesoj
            sompW = sompW + (pg(npj)%pres - pi) * pesoj
         enddo
#ifdef SOLID_BODIES
            sompW = sompW + sompW_vec(npi)
            AppUnity = AppUnity + AppUnity_vec(npi)
#endif
         if (Domain%tipo=="bsph") then
            TetaP1 = input_any_t%TetaP * Med(pg(npi)%imed)%Celerita * dt /     &
                     Domain%h
            pg(npi)%vpres = pi + TetaP1 * sompW / AppUnity
            else
               VIntWdV_FT = zero
               VIntWdV_SO = zero
               VIntWdV_OSB = zero
               VIntWdV_OSP = zero
               press_so = pi
               press_osb = pi
               Ncbs = BoundaryDataPointer(1,npi)
               IntNcbs = BoundaryDataPointer(2,npi)
               ibdt = BoundaryDataPointer(3,npi)
               do icbs=1,IntNcbs
                  ibdp = ibdt + icbs - 1
                  IntLocXY(1:PLANEDIM) =BoundaryDataTab(ibdp)%LocXYZ(1:PLANEDIM)
                  iside = BoundaryDataTab(ibdp)%CloBoNum
                  sidestr = BoundarySide(iside)%stretch
                  strtype = Tratto(sidestr)%tipo
                  do i=1,PLANEDIM
                     sss(i) = BoundarySide(iside)%T(acix(i),1)
                     nnn(i) = BoundarySide(iside)%T(acix(i),3)
                  enddo
                  massforce(1) = sss(1) * Domain%grav(acix(1)) + sss(2) *      &
                                 Domain%grav(acix(2))
                  massforce(2) = nnn(1) * Domain%grav(acix(1)) + nnn(2) *      &
                                 Domain%grav(acix(2))
                  IntWdV = BoundaryDataTab(ibdp)%BoundaryIntegral(2)
                  IntDpWdV(1:2) = BoundaryDataTab(ibdp)%BoundaryIntegral(4:5)
                  if ((strtype=="fixe").or.(strtype=="tapi")) then
                     VIntWdV_FT = VIntWdV_FT + IntWdV
                     sompW = sompW + ro0i * (massforce(1) * IntDpWdV(1) +      &
                             massforce(2) * IntDpWdV(2))
                     elseif (strtype=="sour") then
                        VIntWdV_SO = VIntWdV_SO + IntWdV
                        izonelocal = pg(npi)%izona
                        if (partz(izonelocal)%pressure=="pa") then
                           press_so = partz(izonelocal)%valp
                           elseif (partz(izonelocal)%pressure=="qp") then
                              press_so = ro0i * Domain%grav(3) *               &
                                 (pg(npi)%coord(3) - partz(izonelocal)%valp)
                        endif
                        elseif (strtype=="crit") then
                           VIntWdV_OSP = VIntWdV_OSP + IntWdV
! Implicitly press_osp = pi
                           elseif (strtype=="leve") then
                              VIntWdV_OSB = VIntWdV_OSB + IntWdV
                              izonelocal = pg(npi)%izona
                              vin = BoundarySide(iside)%T(acix(1), 3) *        &
                                    pg(npi)%Vel(acix(1)) +                     &
                                    BoundarySide(iside)%T(acix(2), 3) *        &
                                    pg(npi)%Vel(acix(2))
                              if (vin<zero) then
                                 press_osb = ro0i * Domain%grav(3) *           &
                                    (pg(npi)%coord(3) - partz(izonelocal)%valp)
                                 else
                                    press_osb = pi
                              endif
                  endif
               enddo
! Computing TetaP depending on the time step
               TetaP1 = input_any_t%TetaP * Med(pg(npi)%imed)%Celerita * dt /  &
                        Domain%h
               AppUnity = AppUnity + VIntWdV_FT + VIntWdV_SO + VIntWdV_OSP +   &
                          VIntWdV_OSB
               select case (SmoothVersion)
                  case ("a")
                     smoothpi = pi + TetaP1 * (sompW + (press_so - pi) *       &
                                VIntWdV_SO + (press_osb - pi) * VIntWdV_OSB)   &
                                / AppUnity
                  case ("b")
                     smoothpi = pi + TetaP1 * (sompW + (press_so - pi) *       &
                                VIntWdV_SO + (press_osb - pi) * VIntWdV_OSB +  &
                                (p0i - pi) * (one - AppUnity))
                  case ("c")
                     smoothpi = pi + TetaP1 * (sompW + (press_so - pi) *       &
                                VIntWdV_SO + (press_osb - pi) * VIntWdV_OSB)
                  case default
! no smoothing
                     smoothpi = Pg(npi)%pres       
               endselect
               pg(npi)%vpres = smoothpi
         endif
      endif
   endif
enddo
!$omp end parallel do
! The new density and pressure values are temporarily saved in PartDens and 
! PartPress. In the next cycle, these values are copied in row k.
!$omp parallel do default(none)                                                &
!$omp shared(nag,pg,Domain,Med,indarrayFlu,Array_Flu,input_any_t)              &
!$omp private(npi,ii)
do ii=1,indarrayFlu
   npi = Array_Flu(ii)
! Excluding particles close to the face with conditions "flow", "velo" and 
! "sour"
   if (pg(npi)%koddens==0) then 
      pg(npi)%pres = pg(npi)%vpres
! EoS inverse
      call EoS_barotropic_linear(Med(pg(npi)%imed)%eps,Med(pg(npi)%imed)%den0, &
         Domain%prif,p_in=pg(npi)%vpres,rho_out=pg(npi)%dens)
! Mass update
      if (input_any_t%ALE3) then
         pg(npi)%mass = pg(npi)%dens * pg(npi)%volume
      endif
   endif
enddo
!$omp end parallel do
!------------------------
! Deallocations
!------------------------
#ifdef SOLID_BODIES
   deallocate(sompW_vec)
   deallocate(AppUnity_vec)
#endif
return
end subroutine PressureSmoothing_2D
#endif
