!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: fixed_bed_slope_limited
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Forced deposition or no erosion for particles at least 2h below 
!                  the fixed (as it is defined in the column9 during the same time step: i.e. the maximum 
!                  slope of the fixed bed is 2h/2h. This avoid an eventual too fast propagation of erosion within the column (erosion is a surface phenomena).
!
! Calling routines: Shields
!
! Called subroutines: / 
!
!************************************************************************************

 subroutine fixed_bed_slope_limited(npi,igridi,jgridi,test)

! Assigning modules
 use GLOBAL_MODULE 
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations
 implicit none
 integer(4),intent(in) :: npi,igridi,jgridi
 logical,intent(inout) :: test
 double precision :: Velocity2,fixed_bed_tolerance
 double precision :: aux_vec(3)
 integer(4) :: aux_ID
 double precision :: pretot
!Initializations 
 if (ncord==3) then
     fixed_bed_tolerance = 4.d0*Domain%h
     elseif (ncord==2) then
        fixed_bed_tolerance = 2.d0*Domain%h
 endif    
! Statements
 if (ind_interfaces(igridi,jgridi,4)>0) then
    if (pg(npi)%coord(3)<(pg(ind_interfaces(igridi,jgridi,4))%coord(3)-fixed_bed_tolerance)) then
       if (pg(npi)%state=="flu") then
          pg(npi)%state = "sol"   
          pg(npi)%vel = 0.d0
          pg(npi)%var = 0.d0
          pg(npi)%sigma_prime = 0.0d0  
       endif
       if (pg(npi)%indneighliqsol.ne.0) then
          aux_ID = pg(npi)%indneighliqsol
          elseif (pg(npi)%state=="flu") then 
             aux_ID = pg(npi)%ind_neigh_mob_for_granmob    
             else
                aux_ID = pg(npi)%ind_neigh_mix_bed 
       endif
       if (aux_ID.ne.0) then
          aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
          Velocity2 = dot_product(aux_vec,aux_vec) 
          pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 + Velocity2 * half * Med(pg(aux_ID)%imed)%den0
          pg(npi)%dens = med(pg(npi)%imed)%den0 + (pretot / (Med(pg(npi)%imed)%celerita*Med(pg(npi)%imed)%celerita))   
          else
             if (ind_interfaces(igridi,jgridi,3)>0) then
                pretot = pg(ind_interfaces(igridi,jgridi,3))%pres  + &
                         (pg(ind_interfaces(igridi,jgridi,3))%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(pg(npi)%imed)%den0
                pg(npi)%dens = med(pg(npi)%imed)%den0 + (pretot / (Med(pg(npi)%imed)%celerita*Med(pg(npi)%imed)%celerita))
             endif
       endif 
    test = .false.
    endif 
 endif

 return
 end subroutine fixed_bed_slope_limited
!---split

