!cfile JdWsRn.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : JdWsRn
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
! Module purpose : Module to compute boundary contributions to rodivV
!
! Calling routine: ComputeKernelTable
!
! Called routines: 
!
!************************************************************************************
!
double precision Function JdWsRn(ro, SD, n, kernel)
!
!.. assign modules
  use GLOBAL_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
  integer(4),      intent(IN)      :: SD, n, kernel
  double precision,intent(IN)      :: ro
!
!.. Parameters
  double precision,parameter :: KS2D = 0.454728408833987d0 != 10 / (7*pigreco)
  double precision,parameter :: KA2D = 0.09947192d0        !=5./(16*pigreco)
  double precision,parameter :: KS3D = 0.31830989d0        != 1 / pigreco
  double precision,parameter :: KA3D = 0.07460394d0        !=15./(16*pigreco)
!
!..  0.1875                     = 3/16
!..  0.5625                     = 9/16
!
!.. Local Scalars ..
  double precision :: ro2, ro3, ro4, ro5, duemro
!
!.. Executable Statements ..
!
  ro2 = ro * ro
  ro3 = ro2 * ro
  ro4 = ro2 * ro2
  ro5 = ro4 * ro
!
  JdWsRn = zero
  Select Case (kernel)
    Case (1)              !kernel spline cubica (Monagan)
    
        Select Case (n)
            Case (0)      !n = 0
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -one + (1.5d0 - 0.75d0 * ro) * ro2
                else if (ro >= one .And. ro < two) then
                  duemro = two - ro
                  JdWsRn = -0.25d0 * duemro * duemro * duemro
                end if
            Case (1)      !n = 1
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -0.75d0 + (one - 0.5625d0 * ro) * ro3
                else if (ro >= one .And. ro < two) then
                  JdWsRn = -one + (1.5d0 - one * ro + 0.1875d0 * ro2) * ro2
                end if
            Case (2)      !n = 2
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -0.7d0 + (0.75d0 - 0.45d0 * ro) * ro4
                else if (ro >= one .And. ro < two) then
                  JdWsRn = -0.8d0 + (one - 0.75d0 * ro + 0.15d0 * ro2) * ro3
                end if
            Case (3)      !n = 3
                if (ro >= zero .And. ro < one) then
                  JdWsRn = -0.75d0 + (0.6d0 - 0.375d0 * ro) * ro5
                else if (ro >= one .And. ro < two) then
                  JdWsRn = -0.8d0 + (0.75d0 - 0.6d0 * ro + 0.125d0 * ro2) * ro4
                end if
            case default
                JdWsRn = zero
        End Select
        
        Select Case (SD)
            Case (2)      !geometria 2D
                JdWsRn = JdWsRn * KS2D
            Case (3)      !geometria 3D
                JdWsRn = JdWsRn * KS3D
            case default
                JdWsRn = zero
        End Select
        
    Case (2)              !kernel anticluster cubica (Gallati)
    
        Select Case (n)
            Case (1)      !n = 1
                if (ro < two) then
                  JdWsRn = -4.0d0 + (6.0d0 - 4.0d0 * ro + 0.75d0 * ro2) * ro2
                end if
            Case (2)      !n = 2
                if (ro < two) then
                  JdWsRn = -3.2d0 + (4.0d0 - 3.0d0 * ro + 0.6d0* ro2) * ro3
                end if
            Case (3)      !n = 3
                if (ro < two) then
                  JdWsRn = -3.2d0 + (3.0d0 - 2.4d0 * ro + 0.5d0 * ro2) * ro4
                end if
            case default
                JdWsRn = zero
        End Select
        
        Select Case (SD)
            Case (2)      !geometria 2D
                JdWsRn = JdWsRn * KA2D
            Case (3)      !geometria 3D
                JdWsRn = JdWsRn * KA3D
            case default
                JdWsRn = zero
        End Select

    case default
        
  End Select
!
return
End Function JdWsRn
!---split

