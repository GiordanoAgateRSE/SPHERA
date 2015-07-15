!cfile ReadInputControlSections.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
subroutine ReadInputControlSections ( NumberEntities,Control_Sections, &
                                      ainp,comment,nrighe,ier,ninp,nout )

use GLOBAL_MODULE
use AdM_USER_TYPE

implicit none

integer(4)    :: nrighe,ier, ninp,nout, npts
integer(4),        dimension(20)             :: NumberEntities
type (TySection),  dimension(0:Nsections+1)  :: Control_Sections

character( 1) :: comment
character(80) :: ainp

integer(4)    :: icord, icor2, icor3, icolor, ndiv, ioerr
character(8)  :: label
character(80) :: token
double precision, dimension(3)   :: vp
double precision, dimension(3,2) :: values

character(1), dimension(3) :: CoordLabel = (/ "x", "y", "z" /)

integer(4),    external :: NumberSectionPoints
character(80), external :: lcase, GetToken
logical,       external :: ReadCheck


 call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
 if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL SECTIONS DATA",ninp,nout) ) return

 npts = npoints+npointsl

 do while ( TRIM(lcase(ainp)) /= "##### end control sections #####" )

    NumberEntities(12) = NumberEntities(12) + 1
!   if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION LABEL",ninp,nout) ) return
    label(1:8) = ainp(1:8)

    vp = zero
    values(:,1) = -99999999.
    values(:,2) =  99999999.

   !codice coordinata costante e suo valore
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:1) )
          case ( "x", "y" ,"z" )
             icord = index("xyz",token(1:1))
             token = lcase(GetToken(ainp,2,ioerr))
             if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) vp(icord)
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - CONSTANT COORD. DEFINITION",ninp,nout) ) return

   !codice prima coordinata limitata e valori min. max. (esclusa la coord. costante!)
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:1) )
          case ( "x", "y" ,"z" )
             icor2 = index("xyz",token(1:1))
             if ( icor2 /= icord ) then
                token = lcase(GetToken(ainp,2,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor2,1)
                if ( ioerr == 0 ) token = lcase(GetToken(ainp,3,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor2,2)
             else
                ioerr = -1 !forza errore caso opzione errata
             end if
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - FIRST LIMIT DEFINITION",ninp,nout) ) return

   !codice seconda coordinata limitata e valori min. max. (esclusa la prima coord. e quella costante!)
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:1) )
          case ( "x", "y" ,"z" )
             icor3 = index("xyz",token(1:1))
             if ( icor3 /= icord .AND. icor3 /= icor2 ) then
                token = lcase(GetToken(ainp,2,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor3,1)
                if ( ioerr == 0 ) token = lcase(GetToken(ainp,3,ioerr))
                if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) values(icor3,2)
             else
                ioerr = -1 !forza errore caso opzione errata
             end if
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - SECOND LIMIT DEFINITION",ninp,nout) ) return

   !codice colore ptcl
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    token = lcase(GetToken(ainp,1,ioerr))
    if ( ioerr == 0 ) then
       select case ( token(1:5) )
          case ( "color" )
             token = lcase(GetToken(ainp,2,ioerr))
             if ( ioerr == 0 ) read ( token,*,iostat=ioerr ) icolor
          case default
             ioerr = -1 !forza errore caso opzione NON riconosciuta
       end select
    end if
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONTROL SECTION "//label//" - COLOR INDEX",ninp,nout) ) return
  
   !Calcola il numero di punti da generare
    Ndiv = NumberSectionPoints ( values, CoordLabel(icord) )
    NumberEntities(13) = NumberEntities(13) + Ndiv

    if ( ncord > 0 ) then

       Control_Sections(NumberEntities(12))%Label             = label
       Control_Sections(NumberEntities(12))%Tipo              = CoordLabel(icord)  !//"+"
       Control_Sections(NumberEntities(12))%Constant(:)       = vp(:)
       Control_Sections(NumberEntities(12))%icont(1)          = npts + 1
       Control_Sections(NumberEntities(12))%icont(2)          = npts + Ndiv
       Control_Sections(NumberEntities(12))%XYZRange(:,1)     = -99999999.
       Control_Sections(NumberEntities(12))%XYZRange(:,2)     =  99999999.
       Control_Sections(NumberEntities(12))%XYZRange(icor2,1) = minval(values(icor2,:))
       Control_Sections(NumberEntities(12))%XYZRange(icor2,2) = maxval(values(icor2,:))
       Control_Sections(NumberEntities(12))%XYZRange(icor3,1) = minval(values(icor3,:))
       Control_Sections(NumberEntities(12))%XYZRange(icor3,2) = maxval(values(icor3,:))
       Control_Sections(NumberEntities(12))%TGLsection        = zero
       Control_Sections(NumberEntities(12))%TGLsection(1,1)   = one
       Control_Sections(NumberEntities(12))%TGLsection(2,2)   = one
       Control_Sections(NumberEntities(12))%TGLsection(3,3)   = one
       Control_Sections(NumberEntities(12))%ColorCode         = icolor

       call  CreateSectionPoints ( vp, values, CoordLabel(icord), NumberEntities(12) )
       npts = npts + Ndiv

       if ( nout > 0 ) then
          write (nout,"(1x,a,i3,1x,a)") &
          "Control section   ",NumberEntities(12),"("//Control_Sections(NumberEntities(12))%label//")"
          write (nout,"(1x,a,1x,a,f12.4)") "Constant Coordinate   ", &
          Control_Sections(NumberEntities(12))%Tipo//"=",Control_Sections(NumberEntities(12))%Constant(icord)
          write (nout,"(1x,a,2f12.4)") &
          CoordLabel(icor2)//" Coordinate Limits  ",Control_Sections(NumberEntities(12))%XYZRange(icor2,:)
          write (nout,"(1x,a,2f12.4)") &
          CoordLabel(icor3)//" Coordinate Limits  ",Control_Sections(NumberEntities(12))%XYZRange(icor3,:)
          write (nout,"(1x,a,i12)") &
          "First Point:      ",Control_Sections(NumberEntities(12))%icont(1)
          write (nout,"(1x,a,i12)") &
          "Last  Point:      ",Control_Sections(NumberEntities(12))%icont(2)
       end if

    end if

    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL SECTIONS DATA",ninp,nout) ) return

 end do

return
end subroutine ReadInputControlSections
!---split

