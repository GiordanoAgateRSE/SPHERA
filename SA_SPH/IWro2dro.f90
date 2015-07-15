!cfile IWro2dro.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
double precision Function IWro2dro(ro)

!Computes the definite integral
!2
!S W*(ro')ro2' dro'
!ro
!
!.. assign modules
use AdM_USER_TYPE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Parameters ..
double precision,parameter :: a1 = 0.333333333333333d0 != 1 / 3
double precision,parameter :: a2 = 0.3d0               != 3 / 10
double precision,parameter :: a3 = 0.125d0             != 1 / 8
double precision,parameter :: a4 = 0.25d0
double precision,parameter :: b1 = 0.0625d0             != 1 / 16
double precision,parameter :: b2 = 0.025d0              != 1 / 40
double precision,parameter :: b3 = 4.16666666666667d-03 != 1 / 240
!
!.. Formal Arguments ..
double precision,intent(IN) :: ro
!
!.. Local Scalars ..
double precision :: ro2,ro3,duemro,duemro2,duemro4,IWro2dro1,IWro2dro2
!
!.. Executable Statements ..
!
if (ro >= zero .And. ro < one) then
!    0.333333333333333 = 1 / 3
!    0.3               = 3 / 10
!    0.125             = 1 / 8
  ro2 = ro * ro
  ro3 = ro2 * ro
  IWro2dro1 = (a1 - a2 * ro2 + a3 * ro3) * ro3
  IWro2dro = KERNELCONST3D * (a4 - IWro2dro1)
else if (ro >= one .And. ro < two) then
!    0.0625               = 1 / 16
!    0.025                = 1 / 40
!    4.16666666666667d-03 = 1 / 240
  ro2 = ro * ro
  duemro = two - ro
  duemro2 = duemro * duemro
  duemro4 = duemro2 * duemro2
  IWro2dro2 = -(b1 * ro2 + b2 * duemro * ro + b3 * duemro2) * duemro4
  IWro2dro = KERNELCONST3D * (-IWro2dro2)
Else
  IWro2dro = zero
end if
!
return
End Function IWro2dro
!---split

