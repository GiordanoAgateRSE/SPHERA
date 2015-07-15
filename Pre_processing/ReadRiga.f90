!cfile ReadRiga.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadRiga ( ainp,comment,nrighe,ier,ninp )
implicit none
!
integer(4)   :: ier, ninp
character(1) :: comment
character(*) :: ainp
!
integer(4)   :: ioerr, n, l
integer(4)   :: nrighe
!
  ioerr = 0
!
  READ_LOOP: do while ( ioerr == 0 )
!
    read  ( ninp, "(a)", iostat=ioerr ) ainp
    nrighe = nrighe + 1
!
    if ( ioerr == 0 .AND. trim(ainp) /= "" ) then
!
!.. sostituisce i tabulatori con blank
      if ( ainp(1:1) /= comment ) then
        do n = 1,len(trim(ainp))
          if ( iachar(ainp(n:n)) == 9 ) ainp(n:n) = " "
        end do
        l = index(ainp,comment)
!
        if ( l > 0 ) then
          do n = l,len(trim(ainp))
            ainp(n:n)=" "
          end do
        end if
        exit READ_LOOP 
      end if
!
    end if
!
  end do  READ_LOOP
!
  ier = ioerr
!
return
end subroutine ReadRiga
!---split

