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
! Program unit: OrdGrid1              
! Description: Ordering the numerical elements on the background positioning 
!              grid.       
!-------------------------------------------------------------------------------
subroutine OrdGrid1(nout)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: nout
integer(4) :: npi,ncel,i
integer(4),dimension(Grid%nmax) :: numpartincelgiaaposto
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
numpartincelgiaaposto = 0
Icont = 0
!------------------------
! Statements
!------------------------
! Fluid particles: start
! 1st loop: to find the particle cell and to count the number of particles
! in every cell 
npi = 0
do while ((nag>0).AND.(npi<nag))
   npi = npi + 1
   ncel = ParticleCellNumber(pg(npi)%coord)
   if (pg(npi)%cella<=0) then
      pg(npi) = pg(nag)
      pg(nag) = PgZero
      if (Domain%RKscheme>1) ts0_pg(nag) = ts_pgZero
      nag = nag - 1
      npi = npi - 1
      elseif ((ncel<=0).or.(ncel>Grid%nmax)) then
! The particle is out of grid and is removed from the particle array
         write(nout,'(a,i7,a,i7,3x,3e15.8,a,i7)') "ORDGRID1 particle #",npi,   &
            "   cell:",ncel,pg(npi)%coord(:),                                  &
            "  Total number of particles is = ",nag-1
         EpOrdGrid(pg(npi)%imed) = EpOrdGrid(pg(npi)%imed) + 1
         pg(npi) = pg(nag)
         pg(nag) = PgZero
         if (Domain%RKscheme>1) ts0_pg(nag) = ts_pgZero
         nag = nag - 1
         npi = npi - 1
         else
! To count the particles in each cell and the total number of particles 
! (location nmax+1)
            Icont(ncel) = Icont(ncel) + 1
            pg(npi)%cella = ncel
            Icont(Grid%nmax+1) = Icont(Grid%nmax+1) + 1
   endif
enddo
! 2nd loop: to define the pointer to the first particle in the cell
do i=Grid%nmax,1,-1
   Icont(i) = Icont(i+1) - Icont(i)
enddo
! 3rd loop: counter increased of 1 
do i=1,Grid%nmax+1
   Icont(i) = Icont(i) + 1
enddo
! 4th loop: particle ordering in "NPartOrd"
do npi=1,nag
    ncel = pg(npi)%cella
    if (ncel==0) cycle
    NPartOrd (Icont(ncel) + numpartincelgiaaposto(ncel)) = npi
    numpartincelgiaaposto(ncel) = numpartincelgiaaposto(ncel) + 1
enddo
! Fluid particles: end
! Semi-particles/wall elements (DB-SPH): start
if ((DBSPH%n_w>0).and.((on_going_time_step == it_start).or.                    &
   (Domain%body_part_reorder==1))) then
   Icont_w = 0
   numpartincelgiaaposto = 0
! 1st loop: to find the particle cell and to count the number of particles
! in every cell 
   npi = 0
   do while (npi<(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet))
      npi = npi + 1
      pg_w(npi)%cella = ParticleCellNumber(pg_w(npi)%coord)
! To count the particles in each cell and the total number of particles 
! (location nmax+1)
      Icont_w(pg_w(npi)%cella) = Icont_w(pg_w(npi)%cella) + 1
      Icont_w(Grid%nmax+1) = Icont_w(Grid%nmax+1) + 1
   enddo
! 2nd loop: to define the pointer to the first particle in the cell
   do i=Grid%nmax,1,-1
      Icont_w(i) = Icont_w(i+1) - Icont_w(i)
   enddo
! 3rd loop: counter increased of 1 
   do i=1,Grid%nmax+1
      Icont_w(i) = Icont_w(i) + 1
   enddo
! 4th loop: particle ordering in "NPartOrd_w"
   do npi=1,(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)
NPartOrd_w(Icont_w(pg_w(npi)%cella)+numpartincelgiaaposto(pg_w(npi)%cella))    &
      = npi
      numpartincelgiaaposto(pg_w(npi)%cella) =                                 &
         numpartincelgiaaposto(pg_w(npi)%cella) + 1
   enddo
endif
! Semi-particles/wall elements (DB-SPH): end
! Body particles (Body Transport): start
if (n_bodies>0) then
   Icont_bp = 0
   numpartincelgiaaposto = 0
! 1st loop: to find the particle cell and to count the number of particles
! in every cell
   npi = 0
   do while (npi<n_body_part)
      npi = npi + 1
      bp_arr(npi)%cell = ParticleCellNumber(bp_arr(npi)%pos)
! To count the particles in each cell and the total number of particles 
! (location nmax+1)
      Icont_bp(bp_arr(npi)%cell) = Icont_bp(bp_arr(npi)%cell) + 1
      Icont_bp(Grid%nmax+1) = Icont_bp(Grid%nmax+1) + 1
   enddo
! 2nd loop: to define the pointer to the first particle in the cell
   do i=Grid%nmax,1,-1
      Icont_bp(i) = Icont_bp(i+1) - Icont_bp(i)
   enddo
! 3rd loop: counter increased of 1 
   do i=1,Grid%nmax+1
      Icont_bp(i) = Icont_bp(i) + 1
   enddo
! 4th loop: particle ordering in "NPartOrd_bp"
   do npi=1,n_body_part
NPartOrd_bp(Icont_bp(bp_arr(npi)%cell)+numpartincelgiaaposto(bp_arr(npi)%cell))&
      = npi
      numpartincelgiaaposto(bp_arr(npi)%cell) =                                &
         numpartincelgiaaposto(bp_arr(npi)%cell) + 1
   enddo
endif
! Body particles (Body Transport): end
!------------------------
! Deallocations
!------------------------
return
end subroutine OrdGrid1

