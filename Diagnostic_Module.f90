!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : diagnostic_module
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Added to release 3.1.0
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : provides the global interfaces to the diagnostic routine
!
!************************************************************************************
!
  module diagnostic_module
!
    interface 
      subroutine diagnostic (arg1,arg2,arg3)
        use GLOBAL_MODULE
        integer(4), intent(in)                          :: arg1
        integer(4), intent(in), optional                :: arg2
        character(LEN=lencard), intent(in), optional    :: arg3
      end subroutine 
    end interface 
  end module diagnostic_module
