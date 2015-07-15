!cfile ReadInputControlPoints.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputControlPoints ( NumberEntities,Control_Points, &
                                    ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4),        dimension(20)        :: NumberEntities
type (TyCtlPoint), dimension(NPointst)  :: Control_Points

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,icord,ioerr
double precision, dimension(3) :: values1

character(80), external :: lcase
logical,       external :: ReadCheck


 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL POINTS DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end control points #####" )

    NumberEntities(4) = NumberEntities(4) + 1
    read ( ainp,*,iostat=ioerr ) values1(1:ncord)
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL POINT COORDINATES",ninp,nout) ) return

    if ( ncord > 0 ) then
       Control_Points(NumberEntities(4))%coord(1:3)  = zero
       Control_Points(NumberEntities(4))%dist        = zero
       do n = 1, ncord
          icord = icoordp(n,ncord-1)
          Control_Points(NumberEntities(4))%coord(icord) = values1(n)
       end do
       if ( nout > 0 ) then
          i = NumberEntities(4)
          write (nout,"(1x,a,i5,3(3x,a,1p,e12.4))") &
          "Control point ",i,(xyzlabel(icoordp(n,ncord-1))//" = ",Control_Points(i)%coord(icoordp(n,ncord-1)),n=1,ncord)
       end if
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL POINTS DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputControlPoints
!---split

