!cfile GetToken.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
character(80) function GetToken(ainp,itok,ioerr)
implicit none
!
character(*) :: ainp
integer(4)   :: itok, ioerr
integer(4)   :: n
!
integer(4)   :: number_token
integer(4), dimension(2,20) :: index_token
logical      :: blank 
!
 number_token = 0
 index_token  = 0
 blank        = .TRUE.

 do n = 1, len_trim(ainp)
 
   !if ( ainp(n:n) == "!" ) exit  !Commento dopo i dati
    if ( blank .AND. (ainp(n:n) /= " ") ) then
       number_token = number_token + 1
       index_token(1,number_token) = n
       index_token(2,number_token) = n
       blank = .FALSE.
    else if ( .NOT.blank .AND. (ainp(n:n) /= " ") ) then
       index_token(2,number_token) = n
    else if ( ainp(n:n) == " " ) then 
       blank = .TRUE.
    end if

 end do

 if ( itok <= number_token ) then
    ioerr = 0
    GetToken = ainp( index_token(1,itok):index_token(2,itok) )
 else
    ioerr = itok
    GetToken = ""
 end if
!
return
end function GetToken
!---split

