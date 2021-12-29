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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: wavy_inlet
! Description: To provide a very slightly wavy flow at the inlet section. Each
!              particle layer is staggered by 0.5dx with respect to the previous
!              and the following ones, which are instead aligned each other.
!              This numerical feature reduces the SPH truncation errors at the
!              DB-SPH inlet sections. A white noise is also added. (Amicarelli
!              et al., 2013, IJNME).                
!-------------------------------------------------------------------------------
subroutine wavy_inlet(i_inlet)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_interface_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_inlet
integer(4) :: npart1
#ifdef SPACE_3D
integer(4) :: ii,npart2
#endif
double precision :: rnd
#ifdef SPACE_3D
double precision :: Length1,Length2
#endif
double precision,dimension(3) :: cos_dir_1
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
npart1 = 0
#ifdef SPACE_3D
Length1 = 0.d0
Length2 = 0.d0
npart2 = 0
#endif
!------------------------
! Statements
!------------------------
#ifdef SPACE_2D
   npart1 = NumPartperLine(i_inlet)
   cos_dir_1(:) = BoundarySide(i_inlet)%T(:,1)
#elif defined SPACE_3D
! Length1 and Length 2 are the length scales of the inlet section: they are 
! computed as the distance between the first and the last inlet vertices 
! and the third and the last inlet vertices, respectively.
! Particles are aligned with Plast-P1 and Plast-P3, where P1 is the first 
! boundary vertex, ... and Plast the last boundary vertex. 
! In case of a triangular inlet, we have particles aligned with one direction: 
! P3-P1.
! In case of a quadrilateral inlet, we have particles distributed along two 
! directions: P4-P1 and P4-P3.
      do ii=1,3
         Length1 = Length1 + (BoundaryFace(i_inlet)%Node(1)%GX(ii) -           &
            BoundaryFace(i_inlet)%Node(BoundaryFace(i_inlet)%nodes)%GX(ii)) ** 2
         Length2 = Length2 + (BoundaryFace(i_inlet)%Node(3)%GX(ii) -           &
            BoundaryFace(i_inlet)%Node(BoundaryFace(i_inlet)%nodes)%GX(ii)) ** 2
      enddo
      Length1 = dsqrt(Length1)
      Length2 = dsqrt(Length2)
      npart1 = int(Length1 / Domain%dx + 1.d-2)
      npart2 = int(Length2 / Domain%dx + 1.d-2)
      npart1 = npart1 * npart2 
      cos_dir_1(:) = BoundaryFace(i_inlet)%T(:,1)
#endif
! ID particle (=nag) indicates the last generated particle (of the on-going 
! inlet section)
select case (mod(itime_jet,4))
   case (1,3)  
      pg(nag)%coord(:) = pg(nag)%coord(:) - 0.25d0 * Domain%dx * cos_dir_1(:)
   case (2,0) 
      pg(nag)%coord(:) = pg(nag)%coord(:) + 0.25d0 * Domain%dx * cos_dir_1(:)
endselect
call random_number(rnd)
pg(nag)%coord(:) = pg(nag)%coord(:) + (two * rnd - one) * 0.1d0 * Domain%dx    &
                   * cos_dir_1(:)
!------------------------
! Deallocations
!------------------------
return
end subroutine wavy_inlet
