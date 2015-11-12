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
! Program unit: SetParticles     
! Description: Particle coordinates (initial conditions).                
!----------------------------------------------------------------------------------------------------------------------------------

subroutine SetParticles(Nt,Nz,mate,Xmin,npps,NumParticles,IsopraS)
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
integer(4),intent(IN):: Nt, Nz, mate
integer(4),intent(IN), dimension(SPACEDIM) :: npps
double precision,intent(IN), dimension(SPACEDIM) :: Xmin
integer(4),intent(INOUT) :: NumParticles, IsopraS
logical :: particellainterna
integer(4) :: i,j,k,iaux,test,Nz_aux,nag_aux
double precision :: aux1,aux2,aux3,rnd,tstop
double precision,dimension(SPACEDIM) :: PX
character(len=lencard) :: nomsub = "SetParticles"
logical,external :: IsParticleInternal3D,IsParticleInternal2D
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
   call stoptime (partz(Nz),tstop)
! To compute velocity for particles of type "law"
   call vellaw (partz(Nz)%vlaw,Partz(Nz)%vel,Partz(Nz)%npointv)
end if
if (Domain%tipo=="bsph") then
   if (ncord==3) then
      aux1 = + 0.25d0 * Domain%dd
      else
         aux1 = - 0.25d0 * Domain%dd
   endif
   aux2 = - 0.25d0 * Domain%dd 
   aux3 = - 0.25d0 * Domain%dd 
   iaux = 0
   else
      iaux = 0
      aux1 = - Domain%dd * half
      aux2 = - Domain%dd * half
      aux3 = - Domain%dd * half
endif
PX(1) = Xmin(1) + aux1
! In case the zone is declared but is not used.
if (npps(1)<0) return
! Loop over the X direction.
do i=1,(npps(1)-iaux)
   PX(1) = PX(1) + Domain%dd
   PX(2) = Xmin(2) + aux2
! Loop over the Y direction.
   if (ncord==2) iaux = 0
   do j=1,(npps(2)-iaux)   
      PX(2) = PX(2) + Domain%dd
      PX(3) = Xmin(3) + aux3
! Loop over the Z direction.
      do k=1,(npps(3)-iaux)
         PX(3) = PX(3) + Domain%dd
! To check if the particle falls inside the zone.
         if (ncord==2) then
            particellainterna = IsParticleInternal2D(Nt,PX)
            else 
                particellainterna = IsParticleInternal3D(Nt,PX,IsopraS)
         end if
! In case the particle is inside the domain
         if (particellainterna) then
! the zone counter is increased
            NumParticles = NumParticles + 1
! the total particle number is increased
            if (nagpg==0) cycle
            test = 0
            do Nz_aux=1,NPartZone
               if (Partz(Nz_aux)%IC_source_type==2) test = 1
            end do 
            if (test==0) then
               nag = nag + 1 
! To check the storage for the reached number of particles
               if (nag>PARTICLEBUFFER) call diagnostic                         &
                                               (arg1=10,arg2=4,arg3=nomsub)
               nag_aux = nag 
               else
                  nag_aux = NumParticles 
            endif    
! To modify the coordinates, if random
            if (Domain%RandomPos=='r') then
               call random_number(rnd)
               pg(nag_aux)%coord(1) = PX(1) + (two * rnd - one) * 0.1d0 *      &
                                      Domain%dd
               call random_number(rnd)
               pg(nag_aux)%coord(2) = PX(2) + (two * rnd - one) * 0.1d0 *      &
                                      Domain%dd
               call random_number(rnd)
               pg(nag_aux)%coord(3) = PX(3) + (two * rnd - one) * 0.1d0 *      &
                                      Domain%dd
               else
                  pg(nag_aux)%coord = PX
            end if
            pg(nag_aux)%CoordOld = pg(nag_aux)%coord
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
         end if 
      end do
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine SetParticles

