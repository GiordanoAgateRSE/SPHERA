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
! Program unit: SetParticles    
! Description: Particle coordinates (initial conditions).                
!-------------------------------------------------------------------------------
#ifdef SPACE_3D
subroutine SetParticles(Nz,mate,Xmin,npps,NumParticles,IsopraS)
#elif defined SPACE_2D
subroutine SetParticles(Nz,mate,Xmin,npps,NumParticles)
#endif
!------------------------
! Modules
!------------------------
use Static_allocation_module
use I_O_file_module
use Hybrid_allocation_module
use Dynamic_allocation_module
use Memory_I_O_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: Nz,mate
integer(4),intent(in) :: npps(SPACEDIM)
double precision,intent(in) :: Xmin(SPACEDIM)
integer(4),intent(inout) :: NumParticles
#ifdef SPACE_3D
integer(4),intent(inout) :: IsopraS
#endif
logical :: particellainterna
integer(4) :: i,j,k,iaux,test,Nz_aux,nag_aux,pg_size
double precision :: aux1,aux2,aux3,rnd,tstop
double precision,dimension(SPACEDIM) :: PX
character(len=lencard) :: nomsub = "SetParticles"
#ifdef SPACE_3D
logical,external :: IsParticleInternal3D
#elif defined SPACE_2D
logical,external :: IsParticleInternal2D
#endif
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
if (nagpg>0) then
! To compute "time stop" for particles of type "law"
   call stoptime(partz(Nz),tstop)
! To compute velocity for particles of type "law"
   call vellaw(partz(Nz)%vlaw,Partz(Nz)%vel,Partz(Nz)%npointv)
endif
if (Domain%tipo=="bsph") then
#ifdef SPACE_3D
      aux1 = + 0.25d0 * Domain%dx
#elif defined SPACE_2D
         aux1 = - 0.25d0 * Domain%dx
#endif
   aux2 = - 0.25d0 * Domain%dx 
   aux3 = - 0.25d0 * Domain%dx 
   iaux = 0
   else
      iaux = 0
      aux1 = - Domain%dx * half
      aux2 = - Domain%dx * half
      aux3 = - Domain%dx * half
endif
PX(1) = Xmin(1) + aux1
! In case the zone is declared, but not used.
if (npps(1)<0) return
! Loop over the X direction
do i=1,(npps(1)-iaux)
   PX(1) = PX(1) + Domain%dx
   PX(2) = Xmin(2) + aux2
#ifdef SPACE_2D
      iaux = 0
#endif
! Loop over the Y direction
   do j=1,(npps(2)-iaux)   
      PX(2) = PX(2) + Domain%dx
      PX(3) = Xmin(3) + aux3
! Loop over the Z direction
      do k=1,(npps(3)-iaux)
         PX(3) = PX(3) + Domain%dx
! To check if the particle falls inside the zone
#ifdef SPACE_3D
            particellainterna = IsParticleInternal3D(Nz,PX,IsopraS)
#elif defined SPACE_2D
               particellainterna = IsParticleInternal2D(Nz,PX)
#endif
! In case the particle is inside the zone
         if (particellainterna) then
! The zone counter is increased
            NumParticles = NumParticles + 1
! The total particle number is increased
            if (nagpg==0) cycle
            test = 0
            do Nz_aux=1,NPartZone
               if (Partz(Nz_aux)%IC_source_type==2) test = 1
            enddo
            if (test==0) then
               nag = nag + 1
               nag_aux = nag
               else
                  nag_aux = NumParticles
            endif
! Check the storage for the reached number of fluid particles
            pg_size = size(pg)
            if (nag_aux>pg_size) then
               write(uerr,*) "If you are using a reservoir generated from ",   &
                  "topography, you may need to increase the input parameter ", &
                  "COEFNMAXPARTI. The dimension of the 1D array pg is ",pg_size
               call diagnostic(arg1=10,arg2=4,arg3=nomsub)
            endif
! To modify the coordinates, if random
            if (Domain%RandomPos=='r') then
               call random_number(rnd)
               pg(nag_aux)%coord(1) = PX(1) + (two * rnd - one) * 0.1d0 *      &
                                      Domain%dx
               call random_number(rnd)
               pg(nag_aux)%coord(2) = PX(2) + (two * rnd - one) * 0.1d0 *      &
                                      Domain%dx
               call random_number(rnd)
               pg(nag_aux)%coord(3) = PX(3) + (two * rnd - one) * 0.1d0 *      &
                                      Domain%dx
               else
                  pg(nag_aux)%coord = PX
            endif
            pg(nag_aux)%CoordOld = pg(nag_aux)%coord
            pg(nag_aux)%sect_old_pos(:) = pg(nag_aux)%coord(:)
! Setting Particle Parameters
            if (test==0) then
               call SetParticleParameters(nag,Nz,mate)  
               else
                  call SetParticleParameters(NumParticles,Nz,mate)  
            endif
! Nz is the zone ID, PartZ the vector of zones
            if ((Domain%tipo=="bsph").and.                                     &
               (Partz(pg(nag_aux)%izona)%tipo=="sour"))                        &
               call wavy_inlet(Partz(pg(nag_aux)%izona))
         endif
      enddo
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine SetParticles
