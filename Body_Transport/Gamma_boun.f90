!cfile Gamma_boun.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
! Creation: Amicarelli-Agate 13nov12 
!
!************************************************************************************
! Interpolative function use by Monaghan (2004, 2005) for particle-boundary interactions
double precision function Gamma_boun(r,h)
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
double precision :: r,h,q
!
!.. Executable Statements ..
!
 q = r / h
 if ( q <= (2./3.) ) then
    Gamma_boun = 2./3.
 else if (q <= 1.0) then
    Gamma_boun = 2.*q-(3./2.)*q*q 
 else if (q <=2.0) then
    Gamma_boun = 0.5 * ((2.-q)**2)
 else
    Gamma_boun = 0.
 end if
 
 Gamma_boun = abs(Gamma_boun)

return
end function Gamma_boun 
!---split

