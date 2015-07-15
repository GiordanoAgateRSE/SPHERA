!cfile ReadInputDrawOptions.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputDrawOptions ( ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE                            
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp
character(4)  :: steptime

integer(4)       :: ioerr
character(80)    :: token

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DRAW OPTIONS DATA",ninp,nout) ) return

!do while ( lcase(ainp(1:28)) /= "##### end draw options #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end draw options #####" )

    select case (lcase(GetToken(ainp,1,ioerr)))

       case ("vtkconverter")
          token = lcase(GetToken(ainp,(2),ioerr))
!          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREQUENCY OPTION",ninp,nout) ) return
          select case (token)
             case ("any")
                token = lcase(GetToken(ainp,(3),ioerr))
                read ( token,*,iostat=ioerr ) freq_time
                if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion any :",freq_time," seconds."
                val_time  = zero  !freq_time
             case ("at")
                token = lcase(GetToken(ainp,(3),ioerr))
                read ( token,*,iostat=ioerr ) freq_time
                if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion at :",freq_time," second."
                val_time  = freq_time
                freq_time = -freq_time
             case ("all")
                token = lcase(GetToken(ainp,(3),ioerr))
                read ( token,*,iostat=ioerr ) steptime
                if (steptime == 'time') then
                  freq_time = Domain%memo_fr
                  val_time  = zero  !-Domain%memo_fr
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion every :",freq_time," second."
                else if (steptime == 'step') then
                  freq_time = zero
                  val_time  = const_m_9999
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") "VTKConversion all steps."
                else
                  freq_time = Domain%memo_fr
                  val_time  = zero  !-Domain%memo_fr
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion every :",freq_time," second."
                end if
             case default
                  freq_time = Domain%memo_fr
                  val_time  = zero  !-Domain%memo_fr
                  if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a,1pe12.4,a)") "VTKConversion every :",freq_time," second."
          end select
          if ( ncord > 0 .AND. nout > 0 ) write (nout,"(1x,a)") " "
          vtkconv = .TRUE.

       case default

          ier = 4
          return

    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"DRAW OPTIONS DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputDrawOptions
!---split

