!cfile AddElasticBoundaryReaction_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AddElasticBoundaryReaction_2D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to compute the boundary integral IntWdS
!
! Calling routine: Loop_Irre_2D
!
! Called routines: 
!
!************************************************************************************
!
subroutine AddElasticBoundaryReaction_2D (npi, Ncbs, BoundReaction)

! Adds supplementariìy normal boundary reaction to reinforce insufficient
! pressure gradient in case of few neighbouring particles and presence of
! normal component of mass force (gravity)
! The normal reaction is computed with the formula R=(c0^2/d) ln(zi/d) [for zi<d],
! stemming from the compressible reaction of the fluid, where:
! c0^2 = E/ro0 is the square celerity of the fluid;
! zi is the distance of the particle Pi from the boundary face
! d is a reference distance from which the reaction is added

!.. assign modules
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
integer(4),      intent(IN)    :: npi, Ncbs
double precision,intent(INOUT),dimension(1:SPACEDIM) :: BoundReaction
!
!.. Local Parameters ..
double precision,parameter :: ymincoeff = 0.25d0
double precision,parameter :: reafactor = 1.0d0
!
!.. Local Scalars ..
integer(4) :: sd, icbs, iside, nt, ibdt, ibdp, mate
double precision :: xpi, ypi, ypimin, celer02, vin, normreact
!
!.. Executable Statements ..
!
  mate = pg(npi)%imed
  ypimin = ymincoeff * Domain%dd
  celer02 = Med(mate)%eps / Med(mate)%den0
!
  ibdt = BoundaryDataPointer(3,npi)
  do icbs = 1, Ncbs
!
    ibdp = ibdt + icbs - 1
    iside = BoundaryDataTab(ibdp)%CloBoNum
    nt = BoundarySide(iside)%stretch
!
    if (Tratto(nt)%tipo == "fixe" .or. Tratto(nt)%tipo == "tapi") then
!
      xpi = BoundaryDataTab(ibdp)%LocXYZ(1)
      ypi = BoundaryDataTab(ibdp)%LocXYZ(2)
!
      if (ypi < ypimin) then
!
        if (xpi > zero .and. xpi < BoundarySide(iside)%Length) then  !La proiezione normale della particella npi sul piano
                                                                     ! del lato "iside" è interna al lato "iside"
          vin = zero
          do sd = 1, SPACEDIM
            vin = vin + pg(npi)%var(sd) * BoundarySide(iside)%T(sd, 3)
          end do
          if (vin < zero) then
!            normreact = -reafactor * celer02 * DLog(ypi / ypimin) / ypimin
            normreact = -reafactor * celer02 * DLog((Domain%h + ypi - ypimin) / Domain%h) / Domain%h
            do sd = 1, SPACEDIM
              BoundReaction(sd) = BoundReaction(sd) + normreact * BoundarySide(iside)%T(sd, 3)
            end do
          end if
        end if
      end if
    end if
  end do
!
!prova!
!  BoundReaction = floor(BoundReaction * azzeramento) / azzeramento
!  where (dabs(BoundReaction) < arrotondamento) BoundReaction = zero
!prova!
!
return 
end subroutine AddElasticBoundaryReaction_2D
!---split

