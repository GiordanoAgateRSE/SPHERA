!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-; SPHERA has been authored for RSE SpA by 
!    Andrea Amicarelli, Antonio Di Monaco, Sauro Manenti, Elia Bon, Daria Gatti, Giordano Agate, Stefano Falappi, 
!    Barbara Flamini, Roberto Guandalini, David Zuccalà).
! Main numerical developments of SPHERA: 
!    Amicarelli et al. (2015,CAF), Amicarelli et al. (2013,IJNME), Manenti et al. (2012,JHE), Di Monaco et al. (2011,EACFM). 
! Email contact: andrea.amicarelli@rse-web.it

! This file is part of SPHERA.
! SPHERA is free software: you can redistribute it and/or modify
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
! Program unit: GenerateSourceParticles_2D
! Description: To generate new source particles to simulate inlet fluid flow (only in 2D and with one inlet section).
!----------------------------------------------------------------------------------------------------------------------------------

subroutine GenerateSourceParticles_2D
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use I_O_file_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nt,sd,ip,inttimeratio,isi,i_source
double precision :: Time,SourceTime,TimeFrac,DisplFrac,rnd
character(len=lencard) :: nomsub="GenerateSourceParticles_2D"
integer(4),external :: ParticleCellNumber
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
Time = tempo
!------------------------
! Statements
!------------------------
if ( SourceSide == 0 ) return
inttimeratio = Int(Time / RowPeriod)
if (inttimeratio>pinttimeratio) then
   itime_jet = itime_jet + 1
   SourceTime = inttimeratio * RowPeriod 
   TimeFrac = Time - SourceTime 
   continue
   i_source=0
   do isi=1,NumBSides
      if (BoundarySide(isi)%tipo=="sour") then
         SourceSide = isi
         i_source = i_source + 1
         DisplFrac = RowVelocity(i_source) * TimeFrac 
         do ip=1,NumPartperLine(i_source)
            nag = nag + 1
            SpCount(mat) = SpCount(mat) + 1
            if (nag > PARTICLEBUFFER) then
               call diagnostic (arg1=6,arg2=1,arg3=nomsub)
! To insert message on overpassing the limit or to automatically reallocation 
! of the arrays pg, nPartintorno and associated arrays 
            end if    
! It initializes the parameters of the new particle
            Pg(nag) = PgZero
            if (Domain%RKscheme>1) ts0_pg(nag) = ts_pgZero
            nt = BoundarySide(SourceSide)%stretch
            do sd=1,SPACEDIM
               nn(sd) = BoundarySide(SourceSide)%T(sd,3)
               pg(nag)%coord(sd) = PartLine(i_source,ip,sd) - (yfila +         &
                                   DisplFrac) * nn(sd)
               pg(nag)%vel(sd) = Tratto(nt)%NormVelocity * nn(sd) 
               if (tempo<Tratto(nt)%trampa) pg(nag)%vel(sd) = pg(nag)%vel(sd)  &
                                            * tempo / tratto(nt)%trampa
               pg(nag)%var(sd) = pg(nag)%vel(sd) 
            end do
            if (Domain%tipo == "bsph") call wavy_inlet(isi)
            if (Domain%tipo == "bsph") then
               pg(nag)%rhoSPH_new = zero
               pg(nag)%Gamma = 1.
               pg(nag)%uni = zero
               pg(nag)%sigma = zero
               pg(nag)%dShep = zero 
               pg(nag)%FS = 0 
            endif
            pg(nag)%izona = irz
            pg(nag)%mass = ParticleVolume * Med(Mat)%den0
            pg(nag)%imed = mat           
            pg(nag)%visc = Med(mat)%visc
            pg(nag)%mu = Med(mat)%visc * Med(Mat)%den0
            if ((index(Med(mat)%tipo,"liquid")>0).or.(index(Med(mat)%tipo,     &
               "smagorin")>0)) then
               pg(nag)%state  = "flu"
               pg(nag)%VolFra = VFmn
               else if ((index(Med(mat)%tipo,"granular")>0).or.                &
                  (index(Med(mat)%tipo,"general")>0)) then
                  pg(nag)%state = "sol"
                  pg(nag)%VolFra = VFmx
                  else if (index(Med(mat)%tipo,"gas")>0) then
                     pg(nag)%state = "flu"
                     pg(nag)%VolFra = VFmn
            end if
! Movement/kinematics index
            pg(nag)%vel_type  = partz(irz)%move   
            if (partz(irz)%move/="std") pg(nag)%visc=zero
! Boundary conditions 
            pg(nag)%slip = partz(irz)%slip   
            pg(nag)%cella = ParticleCellNumber(pg(nag)%coord)
! Particle colour 
            call defcolpartzero (irz,partz,pg(nag))
! To initialize "press" and "dens" (density computed from pressure)
            if (partz(irz)%pressure=="pa") then    
! Pression value from input file
               pg(nag)%pres = partz(irz)%valp  
! Free surface level from input file 
               else if (partz(irz)%pressure=="qp") then    
                  pg(nag)%pres = Med(mat)%den0 * Domain%grav(3) *             &
                                 (pg(nag)%coord(3) - partz(irz)%valp)
            end if  
            pg(nag)%dens = Med(mat)%den0 * (one + pg(nag)%pres/Med(mat)%eps)       
         end do
      end if 
   end do 
   pinttimeratio = inttimeratio
end if
!------------------------
! Deallocations
!------------------------
return
end subroutine GenerateSourceParticles_2D
