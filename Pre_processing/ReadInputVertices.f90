!cfile ReadInputVertices.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputVertices ( NumberEntities,Vertice, &
                               ainp,comment,nrighe,ier,prtopt,ninp,nout )

use GLOBAL_MODULE                              
use AdM_USER_TYPE

implicit none

integer(4),      dimension(20)                    :: NumberEntities
double precision,dimension(1:SPACEDIM,NumVertici) :: Vertice

integer(4)    :: nrighe,ier, ninp,nout
logical(4)    :: prtopt
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,icord,ioerr
character(8)  :: label
double precision, dimension(3) :: values1

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end vertices #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,nout) ) return

 if ( ncord > 0 .AND. nout > 0 .AND. prtopt ) then
    write (nout,"(1x,a)") "List of vertices:"
 end if

 do while ( TRIM(lcase(ainp)) /= "##### end vertices #####" )

    select case ( TRIM(Domain%tipo) )
!
!AA406 sub
       case ( "semi","bsph" ) 
!
          read ( ainp,*,iostat=ioerr ) i, values1(1:NumberEntities(1))

! arrotondamento a 1.0d-5
!if ( ncord > 0 ) then
!write (*,'(i10,2f15.10)') i,values1(1:NumberEntities(1))
!values1(1:NumberEntities(1)) = real((nint(values1(1:NumberEntities(1)) * 1.d5)) / 1e5)
!values1(1:NumberEntities(1)) = (anint(values1(1:NumberEntities(1)) * 1.0d5)) / 1.0d5
!write (*,'(i10,2f15.10)') NumberEntities(1),values1(1:NumberEntities(1))
!end if
!

          write (label,"(i8)") i
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTEX n."//label,ninp,nout) ) return
!if ( i > 8 ) i = i - 3420
          NumberEntities(7) = max(i,NumberEntities(7))

          if ( ncord > 0 ) then
             do n = 1, NumberEntities(1)
                icord = icoordp(n,ncord-1)
                if ( NumberEntities(7) == 1 ) then
                   Domain%coord(icord,1) = values1(n)
                   Domain%coord(icord,2) = values1(n)
                end if
                Vertice(icord,i) = values1(n)
                Domain%coord(icord,1) = min(values1(n),Domain%coord(icord,1))
                Domain%coord(icord,2) = max(values1(n),Domain%coord(icord,2))
             end do
          end if

       case default

          if ( nout > 0 ) then
             write (nout,*) "Unknown Domain Type: ",Domain%tipo
          end if
          ier = 2
          return


    end select

    if ( ncord > 0 .AND. nout > 0 .AND. prtopt ) then
       write (nout,"(i6,1p,3(2x,a,e12.4))") &
       i,(xyzlabel(icoordp(n,ncord-1)),Vertice(icoordp(n,ncord-1),i),n=1,ncord)
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"VERTICES DATA",ninp,nout) ) return

 end do

 if ( ncord > 0 .AND. nout > 0 ) then

!   write (nout,*)
!   write (nout,"(1x,a)") "List of vertices:"
!   do n = 1, NumberEntities(7)
!      write (nout,"(i6,1p,3(2x,a,e12.4))") &
!      n,(xyzlabel(icoordp(i,ncord-1)),Vertice(icoordp(i,ncord-1),n),i=1,ncord)
!   end do

    do n = 1, NumberEntities(1)
       icord = icoordp(n,ncord-1)
       write (nout,"(1x,a,a,1p,e12.4)") xyzlabel(icord)," coordinate min. ",Domain%coord(icord,1)
       write (nout,"(1x,a,a,1p,e12.4)") xyzlabel(icord)," coordinate max. ",Domain%coord(icord,2)
    end do
    write (nout,"(1x,a)") " "

 end if

return
end subroutine ReadInputVertices
!---split

