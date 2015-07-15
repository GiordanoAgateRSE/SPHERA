!cfile ltrim.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
character(10) function ltrim(txt)
implicit none
character(*) :: txt

integer(4)   :: i,l,n

 l = len_trim(txt)
 do n = 1,l
    if ( txt(n:n) /= " " ) then
       txt(1:l-n+1) = txt(n:l)
       do i = l-n+2, l
          txt(i:i) = " "
       end do
       exit
    end if
 end do
 ltrim = trim(txt)

return
end function ltrim
!---split

