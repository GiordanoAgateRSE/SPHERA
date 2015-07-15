!cfile ReadInputRestart.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputRestart ( ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                             
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout
logical       :: restartOK
character( 1) :: comment
character(80) :: ainp

integer(4)    :: ioerr
character(80) :: token

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
!  if (restart) then
!    do while ( TRIM(lcase(ainp)) /= "##### end restart #####" )
!      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,nout) ) return
!    end do
!    return
!  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end restart #####" )

    select case ( TRIM(lcase(GetToken(ainp,1,ioerr))) )

       case ("step")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA STEP value",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) Domain%istart
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA STEP value",ninp,nout) ) return
          if ( ncord > 0 .AND. nout > 0 ) then
             write (nout,"(1x,a,i12)") "Restart from step: ",Domain%istart
             if ( Domain%istart < 0 ) write (nout,"(1x,a)") "Negative restart step!"
            !resta attiva solo l'ultima opzione letta
             if ( Domain%start > zero ) then
                write (nout,"(1x,a,f20.12,a)") "Restart from time: ",Domain%start," option ignored!"
                Domain%start = zero
             end if
          end if

       case ("time")
          token = lcase(GetToken(ainp,2,ioerr))
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA TIME value",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) Domain%start
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA TIME value",ninp,nout) ) return
          if ( ncord > 0 .AND. nout > 0 ) then
             write (nout,"(1x,a,f20.12)") "Restart from time: ",Domain%start
             if ( Domain%start < zero ) write (nout,"(1x,a)") "Negative restart time!"
            !resta attiva solo l'ultima opzione letta
             if ( Domain%istart > 0 ) then
                write (nout,"(1x,a,i12,a)") "Restart from step: ",Domain%istart," option ignored!"
                Domain%istart = 0 
             end if
          end if

       case default
          Domain%file = ainp
          if ( ncord > 0 .AND. nout > 0 ) then
             inquire ( file = Domain%file, exist = restartOK )
             if ( restartOK ) then
                write (nout,"(1x,3a)") "Restart file: ",trim(Domain%file)
             else
                write (nout,"(1x,3a)") "Restart file: ",trim(Domain%file)," not found!"
             end if
          end if

    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"RESTART DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputRestart
!---split

