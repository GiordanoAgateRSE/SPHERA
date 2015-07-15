!cfile CreaGrid.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : CreaGrid
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
! Module purpose : Module to create a grid for sorting on the whole domain
!
! Calling routine: Gest_Input
!
! Called routines: diagnostic
!
!************************************************************************************
!
  subroutine CreaGrid
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!  double precision, parameter    :: epsi = 0.001d0
!
!.. local variables
  integer(4) :: ier
  double precision       :: epsi
  character(len=lencard) :: nomsub = "CreaGrid"
!
!.. local arrays
  double precision, dimension(3) :: dextr
!
!.. Executable Statements ..
!
  epsi = 0.01d0 * Domain%dd
!!  epsi = 0.1d0 * Domain%dd
!
!.. set the vertices of the virtual grid on the base of the domain dimensions and the particle smooth,
!.. considering a tolerance epsi on the coordinates for numerical reasons
!
  Grid%extr(:,1) = Domain%coord(:,1) - doubleh - epsi
!
  Grid%extr(:,2) = Domain%coord(:,2) + doubleh + epsi
!
!.. evaluates the dimensions of the grid in all the xyz directions
! 
  dextr(:)       = Grid%extr(:,2) - Grid%extr(:,1)
!
!.. evaluates the number of grid cells in all the directions, assuming that the cell is a cube having 
!.. the side length equal to the double of the smooth length
! 
  Grid%ncd(:)    = NINT(dextr(:) / doubleh)
!
!.. evaluates the real cell dimensions in all the directions, transforming the cubic cell into a real 
!.. hesaedric regular cell
!
  Grid%dcd(:)    = dextr(:) / Grid%ncd(:)
!
!.. in 2D calculation, the number of cells in the Y direction is forced to 1
!
  if ( ncord == 2 ) Grid%ncd(2) = 1
!
!.. evaluates the maximum number of virtual cells in the grid covering the domain parallelepiped
!
  Grid%nmax = Grid%ncd(1) * Grid%ncd(2) * Grid%ncd(3)
!
  write (nout,'(1x,a)') " "
  write (nout,'(1x,a,3i8)') " Number of grid in x, y, z directions : ",Grid%ncd(1),Grid%ncd(2),Grid%ncd(3)
  write (nout,'(1x,a,i10)') " Number of total grid : ",Grid%nmax
  write (nout,'(1x,a)') " "
!
!.. allocazione matrice 2d per calcolo pelolibero (caso erosione)
!
!AA504 sub the fifth element is removed all over the code
  allocate (ind_interfaces(Grid%ncd(1),Grid%ncd(2),4), stat = ier)
  if (ier /= 0) then
!AA504 sub      
    write (nout,'(1x,a,i2)') "    Array ind_interfaces not allocated. Error code: ",ier
    call diagnostic (arg1=4,arg3=nomsub)
  else
!AA504 sub      
    write (nout,'(1x,a)') "    Array ind_interfaces successfully allocated "
  end if
!
  return
  end subroutine CreaGrid
!---split

