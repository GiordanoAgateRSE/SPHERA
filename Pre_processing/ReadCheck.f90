!cfile ReadCheck.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
logical function ReadCheck(IoErr,Ier,Nrighe,ainp,listadati,ninp,nout)
!
!.. implicit declarations
implicit none
!
!.. dummy arguments
integer(4)   :: IoErr,Ier,Nrighe,ninp,nout
character(*) :: ainp, listadati
!
!.. local scalars
!integer(4)   :: MyResult
!character(5) :: Txt
!
  if ( IoErr == 0 ) then
!
    Ier = 0
    ReadCheck = .TRUE.
!
  else
!
    Ier = 4
    ReadCheck = .FALSE.
    write(nout,"(1x,a)")    ">>>>>>>>>>>>>> Warning:"
    write(nout,"(1x,a)")  
    write(nout,"(1x,a,i5)") "Error reading unit:  ",ninp
    write(nout,"(1x,a,a)")  "Expected data:       ",listadati
    write(nout,"(1x,a,i5,a)")  "Last input line read:"//trim(ainp)//"(line number:",Nrighe,")"
!
  end if
!
return
end function ReadCheck
!---split

