!cfile stoptime.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : stoptime
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
! Module purpose : Module for definition of the stop time
!
! Calling routine: SetParticles
!
! Called routines: diagnostic
!
!************************************************************************************
!
  subroutine stoptime ( partzlocal, tstop )
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  type(TyZone),    intent(INOUT)  :: partzlocal
  double precision,intent(INOUT)  :: tstop
!
!.. Local Scalars ..
  integer(4)       :: k,n, icord
  double precision :: tstopc, acc, deltat, spo, dspo, rad
  logical     :: out
  character(len=lencard)  :: nomsub = "stoptime"
!
!.. local arrays
  double precision, dimension(3)   :: dxyz
  double precision, dimension(3,2) :: vlimits, tlimits
!
!.. Executable Statements ..
!
  tstop = max_positive_number
!
!.. checks if there are fixed particles
!
  if ( partzlocal%move == "fix" ) then
!
    do n = 1, ncord
!
       icord = icoordp(n,ncord-1)
       if ( partzlocal%vel(icord) > zero ) then
          tstopc = (Domain%coord(icord,2)-partzlocal%coordMM(icord,2)-one*Domain%dd)/partzlocal%vel(icord)
       else if ( partzlocal%vel(icord) < zero ) then
          tstopc = (Domain%coord(icord,1)-partzlocal%coordMM(icord,1)+one*Domain%dd)/partzlocal%vel(icord)
       else
          tstopc = Domain%tmax
       end if
       tstop = min(tstop,tstopc)
!
     end do
!
!.. checks if there are assigned motion laws for particles
!
   else if ( partzlocal%move == "law" ) then
!
!.. evaluates the paths
!
     vlimits = zero
     do n = 1, ncord
       icord = icoordp(n,ncord-1)
!!       vlimits(icord,1) = Domain%coord(icord,1) + Domain%dd
!!       vlimits(icord,2) = Domain%coord(icord,2) - Domain%dd
       vlimits(icord,1) = Domain%coord(icord,1)
       vlimits(icord,2) = Domain%coord(icord,2)
     end do
!
     tlimits = zero
     dxyz    = zero
     out     = .FALSE.
!
     LAW_ZONE_LOOP: do k = 2, partzlocal%npointv
!
       COORDS_LOOP: do n = 1, ncord
!
         icord = icoordp(n,ncord-1)
!
!.. evaluates the acceleration
!
         deltat = partzlocal%vlaw(0,k) - partzlocal%vlaw(0,k-1)
         acc    = ( partzlocal%vlaw(icord,k) - partzlocal%vlaw(icord,k-1) ) / deltat
!
!.. upgrade the path
!
         dspo = partzlocal%vlaw(icord,k-1) * deltat + acc * deltat * deltat * half
         spo  = dxyz(icord) + dspo
!
!.. checks if the minimum limit has been overridden and how much time has been required
!
         if ( (partzlocal%coordMM(icord,1)+spo) < vlimits(icord,1) ) then
!
           out  = .TRUE.
           dspo = vlimits(icord,1) - (partzlocal%coordMM(icord,1)+dxyz(icord))
           if ( acc == zero ) then
             if (partzlocal%vlaw(icord,k-1) == zero) then
               deltat = max_positive_number
             else
               deltat = dspo / partzlocal%vlaw(icord,k-1)
             end if
           else
             rad  = partzlocal%vlaw(icord,k-1)*partzlocal%vlaw(icord,k-1) - 4.0*0.5*acc*dspo
             if ( rad >= zero ) then
               rad = Dsqrt(rad)  
               deltat = ( partzlocal%vlaw(icord,k-1) + rad ) / ( 2*0.5*acc )
             else
               call diagnostic (arg1=10,arg2=88,arg3=nomsub)       
             end if
           end if
         end if
!
!.. add the interval time to the total time
!
         tlimits(icord,1) = tlimits(icord,1) + deltat  
!
!.. check if the maximum limit has been overriden and how much time has been required 
!
         if ( (partzlocal%coordMM(icord,2)+spo) > vlimits(icord,2) ) then 
           out  = .TRUE.
           dspo = vlimits(icord,2) - (partzlocal%coordMM(icord,2)+dxyz(icord))
           if ( acc == zero ) then
             if (partzlocal%vlaw(icord,k-1) == zero) then
               deltat = max_positive_number
             else
               deltat = dspo / partzlocal%vlaw(icord,k-1)
             end if
           else
             rad  = partzlocal%vlaw(icord,k-1)*partzlocal%vlaw(icord,k-1) - 4.0*0.5*acc*dspo
             if ( rad >= zero ) then
               rad = Dsqrt(rad)  
               deltat = ( -partzlocal%vlaw(icord,k-1) + rad ) / ( 2*0.5*acc )
             else
               call diagnostic (arg1=10,arg2=88,arg3=nomsub)       
             end if
           end if
         end if
!
!.. add the interval time to the total time
!
         tlimits(icord,2) = tlimits(icord,2) + deltat  ! somma il tempo dell'intervallo
!
!.. save the evaluated displacement along the path
!
         dxyz(icord) = dxyz(icord) + dspo
!
       end do COORDS_LOOP
!
!.. ends if it is gone out of the domain
!
       if ( out ) exit LAW_ZONE_LOOP

     end do LAW_ZONE_LOOP
!
!.. evaluates the minimum time
!
     do n = 1, ncord
       icord = icoordp(n,ncord-1)
       tstop = min(tstop,tlimits(icord,1),tlimits(icord,2))
     end do
     partzlocal%move = "fix"
!
   else
!
     tstop = Domain%tmax
!
   end if
!
  return
  end subroutine stoptime
!---split

