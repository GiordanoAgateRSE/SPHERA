!cfile lcase.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
character(80) function lcase(ainp)
implicit none

character(*) :: ainp
integer(4)   :: n, ia

 do n = 1,80
    lcase(n:n) = " "
 end do

 do n = 1, len_trim(ainp)
    ia = iachar(ainp(n:n))
!    write(6,*) n,ia,ainp(n:n); flush(6)
    if ( ia >= 65 .and. ia <= 90 ) then
        lcase(n:n) = char(ia+32)
    else
        lcase(n:n) = ainp(n:n)
    end if
 end do

return
end function lcase
!---split

