!cfile IsPointInternal.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
Logical Function IsPointInternal ( fk, csi )

!Checks wheather a point with local normal coordinates csi() is
!internal to the face whose code is fk (=1 triangle, =2 parallelogram)

!.. assign modules
use GLOBAL_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
integer(4) :: i
integer(4) :: fk
!
!.. Local Arrays ..
double precision, dimension(1:SPACEDIM) :: csi
!
!.. Executable statements ..
!
 IsPointInternal = .FALSE.
!
 if ( fk == 1 ) then            !triangle
   do i = 1, 3
     if ( csi(i) < zero ) return
     if ( csi(i) > one ) return
   end do
   IsPointInternal = .TRUE.
 else if ( fk == 2 ) then        !parallelogram
   do i = 1, 2
     if ( csi(i) < zero ) return
     if ( csi(i) > one ) return
   end do
   IsPointInternal = .TRUE.
 end if
!
return
End Function IsPointInternal
!---split

