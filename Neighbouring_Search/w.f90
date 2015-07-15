!cfile w.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!* kernel standard
double precision function w (r,h,coef)
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a1 = 0.666666667d0
!
!.. Local Scalars ..
double precision :: r,h,s,q,dms,coef
!
!.. Executable Statements ..
!
 s = r / h

 if ( s <= 1.0d0 ) then
    q = a1 + s * s * (s * 0.5d0 - 1.0d0)
 else if (s >= 2.0d0) then
    q = 0.0d0
 else   ! if (s>1.0d0 .and. s<2.0d0) then
    dms = 2.0d0 - s
    q = dms * dms * dms / 6.0d0
 end if

! coef = 0.682093d0/(h*h)

 w = q * coef

return
end function w
!---split

