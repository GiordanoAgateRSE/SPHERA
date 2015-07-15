!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: quadratic_equation
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) To solve a quadratic equation compute k_s 
!
! Calling routines: compute_k_BetaGamma
!
! Called subroutines: / 
!
!************************************************************************************

 subroutine quadratic_equation(a,b,c,n_roots,root1,root2)

! Assigning modules
 use GLOBAL_MODULE 
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations
 implicit none
 double precision, intent(in) :: a,b,c
 integer(4),intent(out) :: n_roots
 double precision, intent(out) :: root1,root2
 double precision :: discriminant
 
! Statements
 root1=-999.d0
 root2=-999.d0
 if (a.ne.0.d0) then
    discriminant = b**2 - 4.0d0*a*c
    if (discriminant>0.0d0) then
       n_roots = 2
       root1 = (-b-dsqrt(discriminant))/(2.0d0*a)
       root2 = (-b+dsqrt(discriminant))/(2.0d0*a)
       else
       if (discriminant==0.0d0) then
          n_roots = 1
          root1 = -b/(2.0d0*a)
          else
             n_roots = 0
       endif      
    endif
    else
       if (b.ne.0.d0) then
          n_roots = 1
          root1 = -c/b
       endif   
  endif       
  
 return
 end subroutine quadratic_equation
!---split

