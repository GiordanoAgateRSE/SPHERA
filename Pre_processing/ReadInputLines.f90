!cfile ReadInputLines.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputLines ( NumberEntities,BoundaryVertex,Tratto, &
                            ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4),               dimension(20)            :: NumberEntities
integer(4),               dimension(NumBVertices)  :: BoundaryVertex
type (TyBoundaryStretch), dimension(NumTratti)     :: Tratto

integer(4)    :: nrighe,ier, ninp,nout
character( 1) :: comment
character(80) :: ainp

integer(4)    :: n,i,ioerr
integer(4)    :: i1,index,numv,numv_line,ipointer
character(5)  :: txt
character(80) :: token

integer(4), parameter :: MAXLINENODES = 20
!integer(4), dimension(MAXLINENODES) :: ivalues

character(80), external :: lcase, GetToken
logical,       external :: ReadCheck

!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end lines #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout) ) return
    end do
    return
  end if
!
 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout) ) return

 do while ( TRIM(lcase(ainp)) /= "##### end lines #####" )

    select case ( TRIM(Domain%tipo) )

!
!AA406 sub
       case ( "semi","bsph" ) 
!

         !Lettura dei vertici che definiscono il boundary
          numv = 0

         !Legge l'indice della linea 
          numv_line = 1
          token = GetToken(ainp,numv_line,ioerr)
          if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) index
          NumberEntities(8) = max(NumberEntities(8),index)  !NumTratti

          VERTEX_LOOP: do while ( ioerr == 0 )
             numv_line = numv_line + 1
             token = GetToken(ainp,numv_line,ioerr)
            !esce quando trova inizio commento o EOR (NON ci sono piu' dati sulla linea di input)
             if ( ioerr /= 0 .OR. &
                  trim(token) == "" .OR. &
                  ichar(trim(token(1:1))) == 9 .OR. &
                  token(1:1) == "!" ) exit VERTEX_LOOP
            !controlla caso codice continuazione linea
             if ( token(1:1) == "&" .OR. token(1:1) == ">" ) then
                call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
                if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VERTICES LIST (continue...)",ninp,nout) ) return
                numv_line = 0
                cycle VERTEX_LOOP
             end if
             numv  = numv + 1
             if ( numv > MAXLINENODES ) then
               stop 'ERRORE in ReadInputLines numv > MAXLINENODES'
             end if
             NumberEntities(9) = NumberEntities(9) + 1  !contatore vertici dei boundaries
             read (token,*,iostat=ioerr) i
             write(txt,"(i5)") i
             if ( numv == 1 ) i1 = i
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VERTEX n."//txt,ninp,nout) ) return
             if ( ncord > 0 ) then
                if ( numv == 1 ) ipointer = NumberEntities(9)
                BoundaryVertex(NumberEntities(9)) = i
             end if
          end do VERTEX_LOOP

         !Conta numero di BoundarySide
          NumberEntities(10) = NumberEntities(10) + numv - 1

          if ( ncord > 0 ) then
             Tratto(index)%numvertices = numv
             Tratto(index)%inivertex   = ipointer
          end if

       case default

          if ( nout > 0 ) then
             write (nout,*) "Unknown Domain Type: ",Domain%tipo
          end if
          ier = 2
          return


    end select

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"LINES DATA",ninp,nout) ) return

 end do


 if ( ncord > 0 .AND. nout > 0 ) then

    write (nout,"(1x,a)") "List of lines"
    write (nout,*)
    do n = 1, NumberEntities(8)
       write (nout,"(1x,a,i3,1x,a)") "Line: ",n
       write (nout,"(1x,a,i3,1x,a)") "Number of Vertices:  ",Tratto(n)%numvertices
       write (nout,"(1x,a,i3,1x,a)") "Vertices Pointer:    ",Tratto(n)%inivertex
       write (nout,"(1x,a,i3,1x,a)") "Vertices List"
       write (nout,"(1x,10i5)") BoundaryVertex(Tratto(n)%inivertex:Tratto(n)%inivertex+Tratto(n)%numvertices-1)
       write (nout,"(1x,a)") " "
    end do

 end if

return
end subroutine ReadInputLines
!---split

