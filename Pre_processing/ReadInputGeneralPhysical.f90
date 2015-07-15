!cfile ReadInputGeneralPhysical.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputGeneralPhysical ( NumberEntities,ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                                     
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4),        dimension(20)     :: NumberEntities

integer(4)    :: n,icord,ioerr
double precision       :: prif
double precision, dimension(3) :: values1

character(80), external :: lcase
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end general physical properties #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",ninp,nout) ) return

!do while ( lcase(ainp(1:43)) /= "##### end general physical properties #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end general physical properties #####" )

    read ( ainp,*,iostat=ioerr ) values1(1:NumberEntities(1))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"GRAVITAL ACCELERATION VECTOR",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    read ( ainp,*,iostat=ioerr ) prif
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"REFERENCE PRESSURE",ninp,nout) ) return

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GENERAL PHYSICAL PROPERTIES DATA",ninp,nout) ) return

 end do

 if ( ncord > 0 ) then
    Domain%grav(:) = zero            
    do n = 1, NumberEntities(1)
       icord = icoordp(n,ncord-1)
       Domain%grav(icord) = values1(n)
       if ( nout > 0 ) write (nout,"(1x,a,a,1p,e12.4)") &
       xyzlabel(icord),"gravity acceler. :",Domain%grav(icord)
    end do
    Domain%prif = prif
    if ( nout > 0 ) write (nout,"(1x,a,1p,e12.4)") "P rif:            :",Domain%prif
    if ( nout > 0 ) write (nout,"(1x,a)") " "
 end if

return
end subroutine ReadInputGeneralPhysical
!---split

