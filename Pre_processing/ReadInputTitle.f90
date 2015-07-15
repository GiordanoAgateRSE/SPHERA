!cfile ReadInputTitle.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputTitle ( ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                            

implicit none

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,ioerr

character(80), external :: lcase
logical,       external :: ReadCheck

 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"TITLE DATA",ninp,nout) ) return 
 n = 0

 do while ( TRIM(lcase(ainp)) /= "##### end title #####" )
 
    n = n + 1
    if ( n <= maxtit ) then
       title(n) = ainp
       if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") title(n)
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"TITLE DATA",ninp,nout) ) return 

 end do
 if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") " "

return
end subroutine ReadInputTitle
!---split

