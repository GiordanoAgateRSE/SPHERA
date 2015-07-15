!cfile ReadInputExternalFile.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!AA504 sub
subroutine ReadInputExternalFile (NumberEntities,ainp,comment,nrighe,ier,OnlyTriangle,ninp,nout,ninp2 )

use GLOBAL_MODULE
use AdM_USER_TYPE
!AA504
use ALLOC_Module

implicit none

integer(4),               dimension(20)                    :: NumberEntities
!AA504 rm part

integer(4)    :: nrighe,ier, ninp,nout, ninp2
character( 1) :: comment
character(80) :: ainp
logical       :: OnlyTriangle

integer(4)    :: ioerr

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end geometry file #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout) ) return

 OnlyTriangle = .TRUE.
 
!do while ( lcase(ainp(1:29)) /= "##### end geometry file #####" )
 do while ( TRIM(lcase(ainp)) /= "##### end geometry file #####" )

    open(ninp2,file=trim(ainp),form="formatted",status="old",iostat=ioerr)

    if ( nout > 0 ) then
       if ( ioerr == 0 ) then
          write (nout,"(1x,3a)") "Geometry File: ",trim(ainp)
       else
          write (nout,"(1x,3a)") "Geometry File: ",trim(ainp)," not found!"
          return
       end if
    end if


   !Legge prima riga file
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp2 )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp2,nout) ) return

    SECTION_LOOP: do while ( ioerr == 0 )

       select case ( TRIM(lcase(trim(ainp))) )

       case ( "##### vertices #####" )

          call ReadInputVertices ( NumberEntities,Vertice, &
                                   ainp,comment,&
                                   nrighe,ier, .FALSE.,ninp2,nout )

       case ( "##### lines #####" )

          call ReadInputLines    ( NumberEntities,BoundaryVertex,Tratto, &
                                   ainp,comment, &
                                   nrighe,ier, ninp2,nout )

       case ( "##### faces #####" )

!AA504 sub           
          call ReadInputFaces    ( NumberEntities,ainp,comment,nrighe,ier,.FALSE.,ninp2,nout)

       case default


       end select

       call ReadRiga ( ainp,comment,nrighe,ioerr,ninp2 )

      !se EOF esce, altrimenti controlla errore
       if ( ioerr == -1 ) cycle SECTION_LOOP

       if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout) ) return

    end do  SECTION_LOOP

    close (ninp2)

    if ( nout > 0 ) then
       write (nout,"(1x,3a)") "End Reading Geometry File"
    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"GEOMETRY FILE",ninp,nout) ) return

 end do

return
end subroutine ReadInputExternalFile
!---split

