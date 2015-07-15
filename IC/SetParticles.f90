!AA504: all the subroutine is modified to treat only coordinates and to call SetParticleParameters
!cfile SetParticles.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : SetParticles
!
! Last updating : April 08, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03  Amicarelli        08Apr14        All the subroutine is modified to treat only coordinates 
!                                      and then to call SetParticleParameters (note: all the comments of v5.03 have been removed)
!AA601
! 04  Amicarelli        26Jan15        DBSPH-input (AA601). DBSPH IC for particle distribution and inlet treatment. 
!
!************************************************************************************
! Module purpose : Module for creation and setting the particles uniformly into the
!                  initial perimeter
!
! Calling routine: GeneratePart
!
! Called routines: IsParticleInternal2D
!                  IsParticleInternal3D
!                  stoptime
!                  vellaw
!AA601
!                  wavy_inlet
!
!
!************************************************************************************

  subroutine SetParticles (Nt, Nz, mate, Xmin, npps, NumParticles, IsopraS)

!Creates and sets particles uniformly into the initial perimeter mib

!Using Modules
  use GLOBAL_MODULE
  use FILES_ENTITIES
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE

!Implicit Declarations
  implicit none

!.. Formal Arguments ..
  integer(4),      intent(IN)                      :: Nt, Nz, mate
  double precision,intent(IN), dimension(SPACEDIM) :: Xmin
  integer(4),      intent(IN), dimension(SPACEDIM) :: npps
  integer(4),      intent(INOUT)                   :: NumParticles, IsopraS

!Local variables
!AA504 sub
  integer(4)       :: i,j,k,iaux,test,Nz_aux,nag_aux
  double precision :: aux1,aux2,aux3,rnd,tstop
  logical          :: particellainterna
  character(len=lencard)  :: nomsub = "SetParticles"
  double precision, dimension(SPACEDIM) :: PX

!External routines
  logical, external    :: IsParticleInternal3D
  logical, external    :: IsParticleInternal2D

!Executable Statements

  if ( nagpg > 0 ) then
!calcolo time stop per particelle tipo 'law'
    call stoptime ( partz(Nz), tstop )
!calcolo velocita' per particelle tipo 'law'
    call vellaw   ( partz(Nz)%vlaw,Partz(Nz)%vel,Partz(Nz)%npointv)
  end if

!AA601 sub start  
  if (Domain%tipo=="bsph") then
     if (ncord==3) then
        aux1 = + 0.25d0*Domain%dd
        else
           aux1 = - 0.25d0*Domain%dd
     endif
     aux2 = - 0.25d0*Domain%dd 
     aux3 = - 0.25d0*Domain%dd 
     iaux = 0
     else
        iaux=0
        aux1 = - Domain%dd * half
        aux2 = - Domain%dd * half
        aux3 = - Domain%dd * half
endif
!AA601 sub end
     
  PX(1) = Xmin(1) +aux1

!In case the zone is declared but is not used
  if (npps(1) < 0) return

!Loops on the X direction

 do i = 1, (npps(1)-iaux)
    PX(1) = PX(1) + Domain%dd
    PX(2) = Xmin(2) + aux2

!Loops on the Y direction
    if (ncord==2) iaux=0
    do j = 1, (npps(2)-iaux)   
       PX(2) = PX(2) + Domain%dd
       PX(3) = Xmin(3) + aux3

!Loops on the Z direction
       do k = 1, (npps(3)-iaux)
          PX(3) = PX(3) + Domain%dd

!Checks if the particle falls inside the zone
          if (ncord == 2) then
             particellainterna = IsParticleInternal2D (Nt, PX)
             else 
                particellainterna = IsParticleInternal3D (Nt, PX, IsopraS)
          end if

!In case the particle is inside the domain
          if ( particellainterna ) then

!the zone counter is increased
             NumParticles = NumParticles + 1

!the total particle number is increased
             if ( nagpg == 0 ) cycle
             
!AA504 sub start 
             test = 0
             do Nz_aux=1,NPartZone
                if (Partz(Nz_aux)%IC_source_type==2) test = 1
             end do 
             if (test==0) then
                nag = nag + 1 
!Check the storage for the reached number of particles
                if (nag > PARTICLEBUFFER) call diagnostic (arg1=10,arg2=4,arg3=nomsub)
                nag_aux = nag 
                else
                   nag_aux = NumParticles 
             endif    
!Modify the coordinates, if random
             if (Domain%RandomPos == 'r') then
               call random_number(rnd)
               pg(nag_aux)%coord(1) = PX(1) + (two * rnd - one) * 0.1d0 * Domain%dd
               call random_number(rnd)
               pg(nag_aux)%coord(2) = PX(2) + (two * rnd - one) * 0.1d0 * Domain%dd
               call random_number(rnd)
               pg(nag_aux)%coord(3) = PX(3) + (two * rnd - one) * 0.1d0 * Domain%dd
               else
                 pg(nag_aux)%coord = PX
             end if
             pg(nag_aux)%CoordOld = pg(nag_aux)%coord
!AA504 sub end             

!Setting Particle Parameters
!AA504 sub start 
             if (test==0) then
                call SetParticleParameters(nag,Nz,mate)  

                else
                   call SetParticleParameters(NumParticles,Nz,mate)  
             endif    
!AA504 sub end             
!AA601 sub start
! Nz is the zone ID, PartZ the vector of zones
             if ((Domain%tipo=="bsph").and.(Partz(pg(nag_aux)%izona)%tipo=="sour")) call wavy_inlet(Partz(pg(nag_aux)%izona))
!AA601 sub end 
          end if !particellainterna
          
        end do
     end do
  end do

  return
  end subroutine SetParticles
!---split

