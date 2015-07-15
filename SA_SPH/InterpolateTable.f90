!cfile InterpolateTable.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine InterpolateTable (xval, nicols, icol, ivalue)
!
!Interpolates values in the array "Table()" with "nrows" rows and "ncols" columns
!Independent variable is in column 0 of Table()
!nicols = number of colums of dependent variables to be interpolated
!icol() = listof columns of dependent variables to be interpolated
!ivalue() = list of the "nicols" interpolated values
!
use Global_MODULE
!
implicit none
!
integer(4)       :: nicols, nr0, nc, ic, nr1
double precision :: xval, xval0, csi, deltaX
!
integer(4),      dimension(1:nicols)         :: icol
double precision,dimension(1:nicols)         :: ivalue
!
deltaX = ktdelta
xval0 = kerneltab(0, 0)
nr0 = Int((xval - xval0) / deltaX)
!
if (nr0 <= 0) then
  Do ic = 1, nicols
    nc = icol(ic)
    ivalue(ic) = kerneltab(0, nc)
  end do
else if (nr0 >= ktrows) then
  Do ic = 1, nicols
    nc = icol(ic)
    ivalue(ic) = kerneltab(ktrows, nc)
  end do
Else
  nr1 = nr0 + 1
  csi = (xval - kerneltab(nr0, 0)) / deltaX
  Do ic = 1, nicols
    nc = icol(ic)
    ivalue(ic) = kerneltab(nr0, nc) + csi * (kerneltab(nr1, nc) - kerneltab(nr0, nc))
  end do
end if
!
return
end subroutine InterpolateTable
!---split

