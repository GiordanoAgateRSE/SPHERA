!cfile NormFix.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine NormFix

use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE

implicit none

integer(4)                    :: npi
double precision              :: unity
double precision,dimension(3) :: appo
!
!.. Executable statements
!
 if ( Ncord == 2 ) then
!
!$omp parallel do default(none) private(npi,appo,unity) shared(nag,Pg)
!
    do npi = 1,nag

       if (pg(npi)%cella == 0 .or. pg(npi)%vel_type == "std") cycle
!
       call InterFix(npi,appo,unity)
!
      !Componenti della normale generica 
       pg(npi)%mno = Dsqrt( (appo(1)*appo(1)) + (appo(3)*appo(3)) )
       pg(npi)%zer(1) = appo(1)/(pg(npi)%mno+0.0001d0)
       pg(npi)%zer(2) = zero
       pg(npi)%zer(3) = appo(2)/(pg(npi)%mno+0.0001d0)
!
       pg(npi)%ang    = (ATAN2 (pg(npi)%zer(1),pg(npi)%zer(3)))
!
       if ( abs(pg(npi)%zer(1)) < 0.5 .AND. abs(pg(npi)%zer(3)) < 0.5 ) then
          pg(npi)%zer(1) = zero
          pg(npi)%zer(2) = zero
          pg(npi)%zer(3) = zero
          pg(npi)%mno    = zero
          pg(npi)%ang    = zero
       end if
!
    end do
!
!$omp end parallel do
!
 else if ( Ncord == 3 ) then
!
!$omp parallel do default(none) private(npi,appo,unity) shared(nag,Pg)
!
    do npi = 1,nag

       if (pg(npi)%cella == 0 .or. pg(npi)%vel_type == "std") cycle

       call InterFix(npi,appo,unity)

       pg(npi)%mno = Dsqrt( (appo(1)*appo(1)) +(appo(2)*appo(2)) + (appo(3)*appo(3)) )
       pg(npi)%zer(:) = appo(:)/(pg(npi)%mno+0.0001d0)

!      pg(npi)%ang  =(ATAN2 (pg(npi)%zer(1),pg(npi)%zer(3)))
!      if ( abs(pg(npi)%zer(1)) < 0.5 .AND. abs(pg(npi)%zer(3)) < 0.5 ) then
   !3D    Decidere che verifica fare in 3D
!         pg(npi)%zer(1) = zero
!         pg(npi)%zer(2) = zero
!         pg(npi)%zer(3) = zero
!         pg(npi)%mno    = zero
!         pg(npi)%ang    = zero
!      end if

    end do
!
!$omp end parallel do
!
 end if
!
return
end subroutine NormFix
!---split

