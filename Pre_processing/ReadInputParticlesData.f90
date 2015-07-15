!cfile ReadInputParticlesData.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
  subroutine ReadInputParticlesData ( NumberEntities, &
                                      Medium,icolor,bends,move,slip,npointv,valuev,values3, &
                                      pressu,valp,ainp,comment,nrighe,ier,ninp,nout )
!
  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
  implicit none
!
  integer(4),        dimension(20)          :: NumberEntities
  integer(4)    :: nrighe,ier, ninp,nout
  character( 1) :: comment
!  character( 8) :: label
  character(80) :: ainp
!  character(4)  :: tipo
!
  integer(4)    :: ioerr
  character(80) :: token
!
  integer(4)    :: Medium
  integer(4)    :: i,n, icord
  integer(4)    :: npointv, icolor
!  double precision       :: trampa
  double precision       :: valp
  character(1)  :: bends
  character(1)  :: slip
  character(3)  :: move
  character(2)  :: pressu
  character(6)  :: token_color
  double precision, dimension(3) :: values3
  double precision, dimension(0:3,maxpointsvlaw) :: valuev
!
  character(80), external :: lcase, GetToken
  logical,       external :: ReadCheck

!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  read ( ainp,*,iostat=ioerr ) Medium
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM INDEX",ninp,nout) ) return
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  token = GetToken(ainp,1,ioerr)
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COLOURING TYPE",ninp,nout) ) return
  bends = lcase(token(1:1))
  token = GetToken(ainp,2,ioerr)
!read ( token,*,iostat=ioerr ) icolor
!
!salva colore per successiva inversione codice
  token_color(1:2) = token(5:6)
  token_color(3:4) = token(3:4)
  token_color(5:6) = token(1:2) 
  read ( token_color,'(Z6)',iostat=ioerr) icolor
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COLOUR",ninp,nout) ) return
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  move = GetToken(ainp,1,ioerr)
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL MOVE STATUS TYPE",ninp,nout) ) return
  values3(:) = zero
!
  slip = " "                             
  npointv = 0
!
  select case (lcase(move))

    case ("std")
       npointv = 0
!       do n = 1, NumberEntities(1)
       do n = 1, 3
          token = GetToken(ainp,(n+1),ioerr)
          read ( token,*,iostat=ioerr ) values3(n)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,nout) ) return
       end do
!       if (tipo == 'sour' .or. tipo == 'velo' .or. tipo == 'flow') then
!         token = GetToken(ainp,(3+2),ioerr)
!         read ( token,*,iostat=ioerr ) trampa
!         if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TEMPO di RAMPA",ninp,nout) ) return
!       end if

    case ("fix")
      !NumberEntities(17) = NumberEntities(17) + 1
       npointv = 1
       valuev  = zero  
       do n = 1, NumberEntities(1)
          token = GetToken(ainp,(n+1),ioerr)
          read ( token,*,iostat=ioerr ) values3(n)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL VELOCITY",ninp,nout) ) return
       end do
      !Aderenza/Free Slip
       token = GetToken(ainp,(NumberEntities(1)+2),ioerr)
       if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NO_SLIP/FREE_SLIP",ninp,nout) ) return
       if      ( lcase(token(1:7)) == "no_slip"   ) then
          slip = "n"
       else if ( lcase(token(1:9)) == "free_slip" ) then
          slip = "f"
       else if ( lcase(token(1:9)) == "cont_slip" ) then  ! 20051212 +2
          slip = "c"
          NumberEntities(19) = 1
       else
          slip = "?"
          ioerr= 1
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NO_SLIP/FREE_SLIP/CONT_SLIP",ninp,nout) ) return
       end if

    case ("law")
       ioerr   = 0
       npointv = 0
       valuev  = zero  
       token = GetToken(ainp,2,ioerr)
       read ( token,*,iostat=ioerr ) npointv
       if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POINTS OF VELOCITY LAW",ninp,nout) ) return
      
       slip = "n"

       do i = 1, npointv
          call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW",ninp,nout) ) return
          if ( ioerr /= 0 .OR. ncord <= 0 ) cycle
          do n = 0, NumberEntities(1)
             icord = icoordp(n,ncord-1)
             token = GetToken(ainp,(n+1),ioerr)
             if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VALUES OF VELOCITY LAW",ninp,nout) ) return
             read ( token,*,iostat=ioerr ) valuev(n,i)
          end do
          if ( i == 1 ) values3(1:3) = valuev(1:3,1)
       end do

    case default

       if ( nout > 0 ) write (nout,*) "Unknown option: ",trim(ainp)
       stop

  end select
!
!.. loads the pressure type and values for all the boundaries
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  pressu = GetToken(ainp,1,ioerr)
  if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PRESSURE TYPE",ninp,nout) ) return
!
!.. the pressure value is assigned ("pa" type)
!
  if ( pressu == "pa" ) then  
!    
    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,nout) ) return
    read ( token,*,iostat=ioerr ) valp
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,nout) ) return
!
!.. the pressure is evaluated as hydrostatic with respect a reference piezo line ("qp" type)
!
  else if ( pressu == "qp" ) then 
!     
    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,nout) ) return
    read ( token,*,iostat=ioerr ) valp
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INITIAL PIEZO LINE",ninp,nout) ) return
!
!.. the pressure is evaluated as hydrostatic with respect the maximum level of an assigned medium ("pl" type)
!
  else if ( pressu == "pl" ) then      
    token = lcase(GetToken(ainp,2,ioerr))
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,nout) ) return
    read ( token,*,iostat=ioerr ) valp
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FREE LEVEL LINE",ninp,nout) ) return
  else
    if ( nout > 0 ) write (nout,*) "Unknown option: ",trim(ainp)
    stop
  end if
!
  if ( ncord > 0 ) then
   !controllo bends
    if ( bends == "b" .AND. icolor > 5 ) then
       icolor = 5
       if ( nout > 0 ) write(nout,*) "Maximum number of bends is 5!"
    end if
  end if
!
  return
!
  end subroutine ReadInputParticlesData
!---split

