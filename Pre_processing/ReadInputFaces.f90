!cfile ReadInputFaces.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!AA504 sub
subroutine ReadInputFaces (NumberEntities,ainp,comment,nrighe,ier,prtopt,ninp,nout)

use GLOBAL_MODULE
use AdM_USER_TYPE
!AA504
use ALLOC_Module

implicit none

integer(4),           dimension(20)        :: NumberEntities
!AA504rm line

integer(4)    :: nrighe,ier, ninp,nout
logical(4)    :: prtopt
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,ioerr, stretch
character(8)  :: label
integer(4), dimension(4) :: ivalues

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck
!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end faces #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end faces #####" )

    select case ( TRIM(Domain%tipo) )

!
!AA406 sub
       case ( "semi","bsph" ) 
!

          ivalues = 0
          read ( ainp,*,iostat=ioerr ) i, ivalues, stretch  ! ivalues(4) - vertice 4 - deve essere 0 se triangolo

          write (label,"(i8)") i
          if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACE n."//label,ninp,nout) ) return

          NumberEntities(11) = max(i,NumberEntities(11))

         !conta i quadrilateri da dividere eventualmente solo su opzione
          if ( ivalues(4)> 0 ) NumberEntities(18) = NumberEntities(18) + 1

          if ( ncord > 0 ) then
             if ( BoundaryFace(i)%Node(1)%name == 0 ) then 
                BoundaryFace(i)%Node(1:MAXFACENODES)%name = ivalues(1:MAXFACENODES)
                BoundaryFace(i)%stretch = stretch
             else
                if ( nout > 0 ) then
                   write (nout,*) "Face definition: ",trim(ainp)
                   write (nout,*) "Face already defined: ",i,BoundaryFace(i)%Node(1:MAXFACENODES)%name
                end if
                ier = 104
                return
             end if
          end if

       case default

          if ( nout > 0 ) then
             write (nout,*) "Unknown Domain Type: ",Domain%tipo
          end if
          ier = 2
          return


    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"FACES DATA",ninp,nout) ) return

 end do


 if ( ncord > 0 .AND. nout > 0 .AND. prtopt ) then

    write (nout,*)
    write (nout,"(1x,a)") "List of faces:"
    do n = 1, NumberEntities(11)
!AA504
       write (nout,"(i10,' - ',4i10,' - ',i8)") n,BoundaryFace(n)%Node(1)%name,BoundaryFace(n)%Node(2)%name,BoundaryFace(n)%Node(3)%name,BoundaryFace(n)%Node(4)%name,BoundaryFace(n)%stretch
    end do

 end if

return
end subroutine ReadInputFaces
!---split

