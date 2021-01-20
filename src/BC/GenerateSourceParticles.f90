!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Description: To generate new source particles at the inlet sections (which 
!              are quadrilaterals in case of 3D simulations). 
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
#ifdef SPACE_3D
   inlet_zone = izone
#elif defined SPACE_2D
      inlet_zone = irz
#endif
!------------------------
! Statements
!------------------------
#ifdef SPACE_3D
   if (SourceFace==0) return
#elif defined SPACE_2D
      if (SourceSide==0) return
#endif
inttimeratio = int(Time / RowPeriod)
if (inttimeratio>pinttimeratio) then
   itime_jet = itime_jet + 1
   SourceTime = inttimeratio * RowPeriod
   TimeFrac = Time - SourceTime
   i_source=0
#ifdef SPACE_3D
      boundary_number = NumFacce
#elif defined SPACE_2D
         boundary_number = NumBSides
#endif
   do isi=1,boundary_number
#ifdef SPACE_3D
         SourceFace = isi
         nt = BoundaryFace(SourceFace)%stretch
         boundary_type = tratto(nt)%tipo
#elif defined SPACE_2D
            boundary_type = BoundarySide(isi)%tipo
#endif
      if (boundary_type=="sour") then
#ifdef SPACE_2D
         SourceSide = isi
#endif
         i_source = i_source + 1
         DisplFrac = RowVelocity(i_source) * TimeFrac
#ifdef SPACE_3D
            NumPartperBound = NumPartFace(i_source)
#elif defined SPACE_2D
               NumPartperBound = NumPartperLine(i_source)
#endif
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
#ifdef SPACE_2D
            nt = BoundarySide(SourceSide)%stretch
#endif
            do sd=1,SPACEDIM
#ifdef SPACE_3D
                  nn(sd) = BoundaryFace(SourceFace)%T(sd,3)
                  pg(nag)%coord(sd) = PartLine(i_source,ip,sd) + (zfila +      &
                                      DisplFrac) * nn(sd)
#elif defined SPACE_2D
                     nn(sd) = BoundarySide(SourceSide)%T(sd,3)
                     pg(nag)%coord(sd) = PartLine(i_source,ip,sd) - (yfila +   &
                                         DisplFrac) * nn(sd)
#endif
               pg(nag)%vel(sd) = Tratto(nt)%NormVelocity * nn(sd)
               pg(nag)%var(sd) = pg(nag)%vel(sd)
            enddo
            pg(nag)%sect_old_pos(:) = pg(nag)%coord(:)
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
            pg(nag)%cella = ParticleCellNumber(pg(nag)%coord)
! Particle colour 
            call defcolpartzero(inlet_zone,partz,pg(nag))
! To initialize "press" and "dens" (density computed from pressure)
            if (partz(inlet_zone)%pressure=="pa") then       
! Pressure value from input file
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
