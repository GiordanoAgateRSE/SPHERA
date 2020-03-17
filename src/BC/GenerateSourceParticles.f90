!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
!--------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: GenerateSourceParticles
! Description: To generate new source particles at the inlet section (only 
!              with one inlet section, which is a quadrilateral in 3D). 
!-------------------------------------------------------------------------------
subroutine GenerateSourceParticles
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
integer(4) :: nt,sd,ip,inttimeratio,isi,i_source,boundary_number,NumPartperBound
integer(4) :: inlet_zone
double precision :: Time,SourceTime,TimeFrac,DisplFrac
character(4) :: boundary_type
character(len=lencard) :: nomsub = "GenerateSourceParticles"
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
Time = simulation_time
if (ncord==3) then
   inlet_zone = izone
   else
      inlet_zone = irz
endif
!------------------------
! Statements
!------------------------
if (ncord==3) then
   if (SourceFace==0) return
   else
      if (SourceSide==0) return
endif
inttimeratio = Int(Time / RowPeriod)
if (inttimeratio>pinttimeratio) then
   itime_jet = itime_jet + 1
   SourceTime = inttimeratio * RowPeriod
   TimeFrac = Time - SourceTime
   i_source=0
   if (ncord==3) then
      boundary_number = NumFacce
      else
         boundary_number = NumBSides
   endif
   do isi=1,boundary_number
      SourceFace = isi
      if (ncord==3) then
         nt = BoundaryFace(SourceFace)%stretch
         boundary_type = tratto(nt)%tipo
         else
            boundary_type = BoundarySide(isi)%tipo
      endif
      if (boundary_type=="sour") then
         if (ncord==2) SourceSide = isi
         i_source = i_source + 1
         DisplFrac = RowVelocity(i_source) * TimeFrac
         if (ncord==3) then
            NumPartperBound = NumPartFace(i_source)
            else
               NumPartperBound = NumPartperLine(i_source)
         endif         
         do ip=1,NumPartperBound
! To generate a new particle row 
            nag = nag + 1
            SpCount(mat) = SpCount(mat) + 1
            if (nag>PARTICLEBUFFER) then         
               call diagnostic(arg1=6,arg2=1,arg3=nomsub)
! To insert message on overpassing the limit or to automatically reallocation 
! of the arrays pg, nPartintorno and associated arrays 
            endif
! Initialization of the quantities of the new particle
            pg(nag) = PgZero
            if (Domain%RKscheme>1) ts0_pg(nag) = ts_pgZero
            if (ncord==2) nt = BoundarySide(SourceSide)%stretch
            do sd=1,SPACEDIM
               if (ncord==3) then
                  nn(sd) = BoundaryFace(SourceFace)%T(sd,3)
                  pg(nag)%coord(sd) = PartLine(i_source,ip,sd) + (zfila +      &
                                      DisplFrac) * nn(sd)
                  else
                     nn(sd) = BoundarySide(SourceSide)%T(sd,3)
                     pg(nag)%coord(sd) = PartLine(i_source,ip,sd) - (yfila +   &
                                         DisplFrac) * nn(sd)
               endif
               pg(nag)%vel(sd) = Tratto(nt)%NormVelocity * nn(sd)
               pg(nag)%var(sd) = pg(nag)%vel(sd)
            enddo
            if (Domain%tipo=="bsph") call wavy_inlet(i_source)
            if (Domain%tipo=="bsph") then
               pg(nag)%rhoSPH_new = zero
               pg(nag)%Gamma = 1.
               pg(nag)%uni = zero
               pg(nag)%sigma = zero
               pg(nag)%sigma_same_fluid = zero
               pg(nag)%dShep = zero 
               pg(nag)%FS = 0 
            endif
! Particle zone
            pg(nag)%izona = inlet_zone
            pg(nag)%mass = ParticleVolume * Med(Mat)%den0
            pg(nag)%imed = mat  
            pg(nag)%kin_visc = Med(mat)%kin_visc
            pg(nag)%mu = Med(mat)%kin_visc * Med(Mat)%den0
            if (index(Med(mat)%tipo,"liquid")>0) then
               pg(nag)%state = "flu"
               elseif (index(Med(mat)%tipo,"granular")>0) then
                  pg(nag)%state = "sol"
            endif
! Movement/kinematics index
            pg(nag)%vel_type = partz(inlet_zone)%move        
            if (partz(inlet_zone)%move/="std") pg(nag)%kin_visc = zero
! Boundary conditions 
            pg(nag)%slip = partz(inlet_zone)%slip        
            pg(nag)%cella = ParticleCellNumber(pg(nag)%coord)
! Particle colour 
            call defcolpartzero(inlet_zone,partz,pg(nag))
! To initialize "press" and "dens" (density computed from pressure)
            if (partz(inlet_zone)%pressure=="pa") then       
! Pression value from input file
               pg(nag)%pres = partz(inlet_zone)%valp + Domain%prif
               elseif (partz(inlet_zone)%pressure=="qp") then   
! Free surface level from input file 
                  pg(nag)%pres = Med(mat)%den0 * Domain%grav(3) *              &
                                 (pg(nag)%coord(3) - partz(inlet_zone)%valp) + &
                                 Domain%prif
            endif
            pg(nag)%dens = Med(mat)%den0 * (one + (pg(nag)%pres - Domain%prif) &
                           / Med(mat)%eps)
         enddo
      endif 
   enddo
   pinttimeratio = inttimeratio
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine GenerateSourceParticles
