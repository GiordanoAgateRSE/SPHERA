!cfile ReadInputDomain.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputDomain ( NumberEntities,ainp,comment,nrighe,ier,ninp,nout,nscr )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout,nscr
character( 1) :: comment
character(80) :: ainp

integer(4),dimension(20) :: NumberEntities

integer(4)       :: ioerr
double precision :: dd, trunc
character(80)    :: token

character(80),external :: lcase, GetToken
logical,      external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end domain #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout) ) return

!do while ( lcase(ainp(1:22)) /= "##### end domain #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end domain #####" )


    token = lcase(GetToken(ainp,1,ioerr))
    read ( token,*,iostat=ioerr ) NumberEntities(1)
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN COORDINATES NUMBER",ninp,nout) ) return 
    if ( ncord > 0 .AND. nout > 0 ) then
       write (nout,"(1x,a,i3,1x,a)") "Domain Dimension       : ",ncord,ncordlabel(ncord)
    end if

    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN TYPE",ninp,nout) ) return 

    select case ( token(1:4) )

!
!AA406 sub
       case ( "bsph","semi" ) 
!
         if (.not. DBSPH_Module .and. token(1:4) == "bsph") then
           write (nscr,"(1x,a)") " "
           write (nout,"(1x,a)") " "
           write (nscr,"(1x,a)") " >>WARNING! - The DBSPH module is not available."
           write (nout,"(1x,a)") " >>WARNING! - The DBSPH module is not available."
           ier = 40
           return
         end if
!
         Domain%tipo = token(1:4)
         if ( ncord > 0 .AND. nout > 0 ) then
           write (nout,"(1x,a,1x,a)"   ) "Domain Type            : ",trim(token)
         end if      

       case default
         ier = 3
         return

    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( ioerr == 0 ) read ( ainp,*,iostat=ioerr ) dd, trunc
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DD & TRUNC",ninp,nout) ) return
    Domain%dd    = dd
    Domain%trunc = trunc

    token = lcase(GetToken(ainp,3,ioerr))
!    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RANDOM INITIAL POSITION",ninp,nout) ) return 
    if (token(1:1) == 'r') then
      Domain%RandomPos = 'r'
    else
      Domain%RandomPos = 'n'
    end if

    if ( ncord > 0 .AND. nout > 0 ) then
       write (nout,"(1x,a,1pe12.4)") "Dd                     : ",dd
       write (nout,"(1x,a,1pe12.4)") "Trunc                  : ",trunc
       write (nout,"(1x,a,1x,a)") "Random Initial Position: ",Domain%RandomPos
       write (nout,"(1x,a)")  " "
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DOMAIN DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputDomain
!---split

