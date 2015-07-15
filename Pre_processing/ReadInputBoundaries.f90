!cfile ReadInputBoundaries.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
  subroutine ReadInputBoundaries ( NumberEntities,Partz,Tratto,BoundaryVertex,ainp,comment,  &
                                   nrighe,ier,ninp,nout )

  use GLOBAL_MODULE
  use AdM_USER_TYPE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
  integer(4),               dimension(20)           :: NumberEntities
  type (TyZone),            dimension(NPartZone)    :: Partz
  type (TyBoundaryStretch), dimension(NumTratti)    :: Tratto
  integer(4),               dimension(NumBVertices) :: BoundaryVertex
  integer(4)                                        :: nrighe,ier, ninp,nout
  character( 1)                                     :: comment
  character(80)                                     :: ainp
!
!.. local scalars

!AA504 sub start
  integer(4)        :: n,index,numv,indexi,indexf,Izona,ipointer,Medium,icolor,icord,ioerr,npointv,IC_source_type,Car_top_zone,dx_CartTopog,plan_reservoir_points,nag_aux
  integer(4)        :: i,i1,i2,i_point,ID_first_vertex,ID_last_vertex,dam_zone_ID,dam_zone_n_vertices
  double precision  :: pool_value,shear,velocity,trampa,valp,flowrate,H_res
!AA504 sub end
  character(len=1)  :: pool_plane,bends,slip
  character(len=2)  :: pressu
  character(len=3)  :: move
  character(len=4)  :: tipo
  character(len=6)  :: token_color
  character(len=8)  :: label
  character(len=80) :: token
!
!.. local arrays
  double precision, dimension(3) :: values1,values3
  double precision, dimension(0:3,maxpointsvlaw) :: valuev
!AA504
  double precision :: plan_reservoir_pos(4,2),dam_zone_vertices(4,2)

!.. external assignments
  integer(4),    external :: ptcolorrgb
  character(80), external :: lcase, GetToken
  logical,       external :: ReadCheck
!
!.. Executable statements
!
!.. initializations
!
  npointv = 0
  values3 = zero
  valp = zero
!
!.. in case of restart the cards are not read
!
  if (restart) then
    do while ( TRIM(lcase(ainp)) /= "##### end boundaries #####" )
      call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
      if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARIES DATA",ninp,nout) ) return
    end do
    return
  end if
!
!.. input loading of boundary set starts...
!
  call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
  if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARIES DATA",ninp,nout) ) return
!
!.. loads all the cards of the set
!
  do while ( TRIM(lcase(ainp)) /= "##### end boundaries #####" )
!
!.. assign the label of the boundary condition
!
    label = ainp(1:8)
!
!.. read the boundary index card for the first zone having the same condition
!
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARIES INDEX",ninp,nout) ) return
    token = GetToken(ainp,1,ioerr)
    read ( token,*,iostat=ioerr ) indexi
!
!.. read the boundary index card for the last zone having the same condition
!
    if ( .NOT.ReadCheck (ioerr,ier,nrighe,ainp,"BOUNDARY INDEX",ninp,nout) ) return
    indexf = indexi
    token = GetToken(ainp,2,ioerr)
    if ( token /= "" )   read ( token,*,iostat=ioerr ) indexf
    NumberEntities(8) = max(indexf,NumberEntities(8))
!
!.. read the boundary type
!
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
    if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BOUNDARY TYPE",ninp,nout) ) return
    tipo = lcase(ainp(1:4))
!
!.. loads the data for the different boundary conditions
!
    numv     = 0
    ipointer = 0
    move     = "   "
    Medium   = 0
    icolor  = Z'FF0000'  
    values1 = zero
    values3 = zero
    pool_plane = " "
    pool_value = zero
    shear      = zero
    velocity   = zero
    flowrate   = zero
    trampa     = zero
    pressu     = "  "
    valp       = zero

!AA504 start
    IC_source_type = 0
    Car_top_zone = 0
    dx_CartTopog = 0.
    H_res = 0.
    ID_first_vertex = 0
    ID_last_vertex = 0
    plan_reservoir_points = 0
    nag_aux = 0
    dam_zone_ID = 0
    dam_zone_n_vertices = 0
    plan_reservoir_pos = 0. 
    dam_zone_vertices = 0.
!AA504 end
    
!    
    select case ( tipo )
!
!.. boundary condition "leve", "crit" or "open"
!
      case ( "leve", "crit", "open" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
          token = GetToken(ainp,1,ioerr)
          token_color(1:2) = token(5:6)
          token_color(3:4) = token(3:4)
          token_color(5:6) = token(1:2) 
          read ( token_color,'(Z6)',iostat=ioerr) icolor
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
!.. boundary condition "fixe"
!
     case ( "fixe" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!

        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) shear
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: SHEAR STRESS COEFFICIENT",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
!.. boundary condition "sour"
!
     case ( "sour" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) Medium
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: MEDIUM INDEX",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: NORMAL VELOCITY, TRAMPA ",ninp,nout) ) return
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE, TRAMPA ",ninp,nout) ) return
!
        token = GetToken(ainp,1,ioerr)
!        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: NORMAL VELOCITY",ninp,nout) ) return
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE",ninp,nout) ) return
!        read ( token,*,iostat=ioerr ) velocity
!        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: NORMAL VELOCITY",ninp,nout) ) return
        read ( token,*,iostat=ioerr ) flowrate
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE",ninp,nout) ) return
        token = GetToken(ainp,2,ioerr)
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: TRAMPA",ninp,nout) ) return
        read ( token,*,iostat=ioerr ) trampa
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: TRAMPA",ninp,nout) ) return
        if (trampa /= 0) then
          write (nout,*) ' '
          write (nout,*) 'TRAMPA in SOURCE boundary is not available. TRAMPA is setted to zero; check the VELOCITY boundary.'
          write (nout,*) ' '
          write (*,*) ' '
          write (*,*) 'TRAMPA in SOURCE boundary is not available. TRAMPA is setted to zero;  check the VELOCITY boundary.'
          write (*,*) ' '
          trampa = zero
        end if
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!        token = GetToken(ainp,1,ioerr)
!        read ( token,*,iostat=ioerr ) pressu
!        pressu = pressu(1:len_trim(pressu))
        pressu = GetToken(ainp,1,ioerr)
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE TYPE",ninp,nout) ) return
        if (pressu == "pa") then  
          token = GetToken(ainp,2,ioerr)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE VALUES",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) valp
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE VALUES",ninp,nout) ) return
        else if (pressu == "qp") then
          token = GetToken(ainp,2,ioerr)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PIEZO LINE",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) valp
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PIEZO LINE",ninp,nout) ) return
!        else if (pressu == "pl") then      
!          token = GetToken(ainp,2,ioerr)
!          read ( token,*,iostat=ioerr ) valp
        else
          if ( nout > 0 ) write (nout,*) "Unknown option: ",trim(ainp)," in source boundary."
          stop
        end if
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
        move   = "std"
!
!.. boundary condition "velo"
!
     case ( "velo" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) velocity,trampa
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VELO: NORMAL VELOCITY, TRAMPA",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
        move   = "std"
        pressu = "pa"
        valp = zero
!
!.. boundary condition "flow"
!
     case ( "flow" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) flowrate,trampa
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"VELO: FLOW RATE, TRAMPA",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z6)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,nout) ) return
!
        move   = "std"
        pressu = "pa"
        valp = zero
!
!.. boundary condition "peri"
!
      case ( "peri" )    
!
        NumberEntities(3) = NumberEntities(3) + 1
!
       
        call ReadInputParticlesData ( NumberEntities, Medium,icolor,bends,move,slip,npointv,valuev,values3, &
                                      pressu,valp,ainp,comment,nrighe,ier,ninp,nout )

        if ( ier /= 0 ) return  
        
!AA504 start
        call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
        if (ioerr==0) read (ainp,*,iostat=ioerr) IC_source_type,Car_top_zone
        if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC_source_type, Car_top_zone",ninp,nout) ) return
        if (IC_source_type==2) then
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) dx_CartTopog,H_res
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"dx_CartTopog,H_res",ninp,nout) ) return
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) ID_first_vertex,ID_last_vertex
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_first_vertex,ID_last_vertex",ninp,nout) ) return
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) plan_reservoir_points,nag_aux
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"plan_reservoir_points and nag_aux",ninp,nout) ) return
           do i2=1,plan_reservoir_points
              call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
              if (ioerr==0) read (ainp,*,iostat=ioerr) plan_reservoir_pos(i2,1),plan_reservoir_pos(i2,2)
              if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"plan_reservoir_vertices",ninp,nout) ) return
           end do
           call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
           if (ioerr==0) read (ainp,*,iostat=ioerr) dam_zone_ID,dam_zone_n_vertices
           if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"dam_zone_ID and dam_zone_vertices",ninp,nout) ) return
           if (dam_zone_ID>1) then
              do i2=1,dam_zone_n_vertices
                 call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
                 if (ioerr==0) read (ainp,*,iostat=ioerr) dam_zone_vertices(i2,1),dam_zone_vertices(i2,2)
                 if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"dam zone vertices",ninp,nout) ) return
              end do
           endif
        endif
!AA504 end
        
!
!.. boundary condition "tapi"
!
      case (  "tapi" )    
!
!.. returns an error if the number of vertices is not equal to two in 3D
!     
        if ( numv /= 2 .and. NumberEntities(1) == 3 ) then 
          write (nout,'(a,i15)') "TAPIS boundary type: 2 vertices are requested:",numv
          ier = 103
          return
        end if
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( ioerr == 0 ) read (ainp,*,iostat=ioerr) shear
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: SHEAR STRESS COEFFICIENT",ninp,nout) ) return
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: VELOCITY COMPONENTS",ninp,nout) ) return
        do n = 1, NumberEntities(1)
          icord = icoordp(n,NumberEntities(1)-1)
          token = GetToken(ainp,n,ioerr)
          if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,xyzlabel(icord)//" VELOCITY COMPONENT (TAPIS)",ninp,nout) ) return
          read ( token,*,iostat=ioerr ) values1(icord)
        end do
!
        call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
        token = GetToken(ainp,1,ioerr)
        token_color(1:2) = token(5:6)
        token_color(3:4) = token(3:4)
        token_color(5:6) = token(1:2) 
        read ( token_color,'(Z)',iostat=ioerr) icolor
        if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: RRGGBB COLOR",ninp,nout) ) return
!
!.. boundary condition "pool" : active only for 3D case
!
       case ("pool")
!
         if ( NumberEntities(1) == 3 ) then      
!   
           NumberEntities(3) = NumberEntities(3) + 1
!
           call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
           pool_plane = GetToken(ainp,1,ioerr)
           if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POOL: X/Y/Z/ PLANE LABEL",ninp,nout) ) return
           token = GetToken(ainp,2,ioerr)
           if ( ioerr == 0 ) read (token,*,iostat=ioerr) pool_value
           if ( .NOT.ReadCheck(ioerr,ier,nrighe,ainp,"POOL: PLANE VALUE",ninp,nout) ) return
!
           call ReadInputParticlesData ( NumberEntities, &
                                         Medium,icolor,bends,move,slip,npointv,valuev,values3, &
                                         pressu,valp,ainp,comment,nrighe,ier,ninp,nout )

           if ( ier /= 0 ) return  
!
         end if 
!
       case default

         write (nout,*) "Unrecognised BOUNDARY type: ",tipo
         ier = 101
         return
!
    end  select
!
!.. load the parameter structure with input data
!
    if ( ncord > 0 ) then
!
       Izona = NumberEntities(3)
!
       Partz(Izona)%label    = label
       Partz(Izona)%tipo     = tipo
       Partz(Izona)%Medium   = Medium

!AA504 start
       Partz(Izona)%IC_source_type = IC_source_type
       Partz(Izona)%Car_top_zone = Car_top_zone
       if (IC_source_type == 2) then
          Partz(Izona)%dx_CartTopog = dx_CartTopog
          Partz(Izona)%H_res = H_res
          Partz(Izona)%ID_first_vertex = ID_first_vertex
          Partz(Izona)%ID_last_vertex = ID_last_vertex
          Partz(Izona)%plan_reservoir_points = plan_reservoir_points
          Partz(Izona)%nag_aux = nag_aux
          Partz(Izona)%plan_reservoir_pos = plan_reservoir_pos
          Partz(Izona)%dam_zone_ID = dam_zone_ID
          Partz(Izona)%dam_zone_n_vertices = dam_zone_n_vertices
          if (dam_zone_ID>1) Partz(Izona)%dam_zone_vertices = dam_zone_vertices
       endif
!AA504 end

       Partz(Izona)%icol     = icolor
       Partz(Izona)%bend     = bends
       Partz(Izona)%move     = move
       Partz(Izona)%slip     = slip
       if (npointv /= 0) then
         Partz(Izona)%npointv = npointv
         Partz(Izona)%vlaw(icoordp(0:ncord,ncord-1),1:npointv)= valuev(0:ncord,1:npointv)
       end if
       Partz(Izona)%vel      = zero
       Partz(Izona)%vel(icoordp(1:ncord,ncord-1)) = values3(1:ncord)
       Partz(Izona)%trampa   = trampa
       Partz(Izona)%pressure = pressu
       Partz(Izona)%valp     = valp
       Partz(Izona)%Indix(1) = indexi
       Partz(Izona)%Indix(2) = indexf
       if ( pool_plane == "X" .OR. pool_plane == "x" ) Partz(Izona)%ipool = 1
       if ( pool_plane == "Y" .OR. pool_plane == "y" ) Partz(Izona)%ipool = 2
       if ( pool_plane == "Z" .OR. pool_plane == "z" ) Partz(Izona)%ipool = 3
       Partz(Izona)%pool = pool_value
!
!.. load the constraints 
!
       MULTI_INDEX_LOOP: do index = indexi, indexf       
!
         Tratto(index)%tipo        = tipo
         if (ncord == 3) then
           Tratto(index)%numvertices = numv
           Tratto(index)%inivertex   = ipointer
         end if
         Tratto(index)%ShearCoeff  = shear
         Tratto(index)%Medium      = Medium
         Tratto(index)%velocity    = values1
         Tratto(index)%NormVelocity= velocity
         Tratto(index)%FlowRate    = flowrate
         Tratto(index)%trampa      = trampa
         Tratto(index)%zone        = Izona
         Tratto(index)%ColorCode   = icolor
!
         if (nout > 0 .and. index == indexi ) then
!
             if ( index > 1 ) write (nout,*)
             if ( indexf == indexi ) write (nout,"(1x,a,i5,1x,a)")    "Boundary        : ",indexi
             if ( indexf /= indexi ) write (nout,"(1x,a,i5,1x,a,i5)") "Boundary        : ",indexi,"   to",indexf
             write (nout,"(1x,a,2x,a)")    "Type            : ",Tratto(index)%tipo 
             if ( tipo == "fixe" ) then
                write (nout,"(1x,a,1pe12.4)") "Shear coeff.    : ",Tratto(index)%ShearCoeff
             else if ( tipo == "peri" ) then
                write (nout,"(1x,a,i3,1x,a)") "Medium Index    : ",Tratto(index)%Medium
             else if ( tipo == "pool" ) then
                write (nout,"(1x,a,i3,1x,a)") "Medium Index    : ",Tratto(index)%Medium
             else if ( tipo == "tapi" ) then
                write (nout,"(1x,a,1pe12.4)") "Shear coeff.    : ",Tratto(index)%ShearCoeff
                do n = 1, ncord
                   icord = icoordp(n,ncord-1)
                   write (nout,"(1x,a,a,1pe12.4)") &
                   xyzlabel(icord)," Velocity      : ",Tratto(index)%velocity(n)
                end do
             end if
             if (ncord == 2) then
               write (nout,"(1x,a)") "Vertices List"
               write (nout,"(1x,10i5)") BoundaryVertex(Tratto(index)%inivertex:Tratto(index)%inivertex+Tratto(index)%numvertices-1)
               write (nout,"(1x,a,z6)") "Color           : ",Tratto(index)%colorCode
             end if
!
          end if
!
          select case (tipo)
!
            case ("fixe")
              if ( ncord == 3) then
                Tratto(index)%ColorCode = icolor
                if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
                Tratto(index)%ColorCode = icolor
                write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
              end if
!
            case ("tapi")
              if (ncord == 3) then 
                Tratto(index)%ColorCode = icolor
                if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
                Tratto(index)%ColorCode = icolor
                write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
              else 
                numv = Tratto(index)%numvertices
                if ( numv /= 2  ) then 
                  if ( nout > 0 ) write (nout,'(a,i15)') "TAPIS boundary type: 2 vertices are requested:",numv
                  ier = 103
                  return
                end if
              end if
!
            case ("peri")
              if (ncord == 2) then 
                i1= BoundaryVertex(Tratto(index)%inivertex)
                i = BoundaryVertex(Tratto(index)%inivertex+Tratto(index)%numvertices-1)
                if ( i /= i1 ) then ! errore se primo ed ultimo vertice sono diversi
                  if ( nout > 0 ) write (nout,'(a,2i15)') "PERIMETER boundary type: first and last vertices are different: ",i,i1
                  ier = 102
                  return
                end if

                if ( nout > 0 ) then

                  write (nout,"(1x,a,i3,1x,a)")  "Zone            : ",Izona,Partz(Izona)%label
                  write (nout,"(1x,a,i3)")       "Medium Index    : ",Partz(Izona)%Medium
                  write (nout,"(1x,a,Z6.6)")     "Color           : ",Partz(Izona)%icol
                  write (nout,"(1x,a,2x,a)")     "Bends           : ",Partz(Izona)%bend
                  write (nout,"(1x,a,2x,a)")     "Movement Type   : ",Partz(Izona)%move
                  write (nout,"(1x,a,2x,a)")     "Boundary Cond.  : ",Partz(Izona)%slip
                  if ( Partz(Izona)%move == "law" ) then
                    write (nout,"(1x,a,i3)")      "Velocity Table - Number of Points: ",Partz(Izona)%npointv
                    do i = 1,Partz(Izona)%npointv
                      write (nout,"(a,i3,1p,4(2x,a,e12.4))") &
                            " Point",i,(xyzlabel(icoordp(n,ncord-1)),Partz(Izona)%vlaw(icoordp(n,ncord-1),i),n=0,ncord)
                    end do
                  end if
                  do n = 1, ncord
                    icord = icoordp(n,ncord-1)
                    write (nout,"(1x,a,a,1pe12.4)") &
                    xyzlabel(icord)," velocity       : ",Partz(Izona)%vel(icord) 
                  end do
                  write (nout,"(1x,a,1pe12.4)") "Time Rampa      : ",Partz(Izona)%trampa
                  write (nout,"(1x,a,2x,a)")    "Pressure Type   : ",Partz(Izona)%pressure
                  write (nout,"(1x,a,1pe12.4)") "Pressure Value  : ",Partz(Izona)%valp
!                 write (nout,"(1x,a,1p,e12.4)") "Monaghan Coeff. : ",Med(Tratto(index)%Medium)%alfaMon
!                 do n = 1, ncord
!                   icord = icoordp(n,ncord-1)
!                   write (nout,"(1x,a,a,1pe12.4)") &
!                   xyzlabel(icord)," min        ",Partz(Izona)%coord(icord,1)
!                   write (nout,"(1x,a,a,1pe12.4)") &
!                   xyzlabel(icord)," max        ",Partz(Izona)%coord(icord,2)
!                 end do
          
                end if ! aggiunta       

              else
                Tratto(index)%ColorCode = icolor
                if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode
              end if
!AA504 start              
              write (nout,"(1x,a,i12)")        "IC_source_type  : ",Partz(Izona)%IC_source_type
              write (nout,"(1x,a,i12)")        "Car_top_zone    : ",Partz(Izona)%Car_top_zone
              if (IC_source_type == 2) then
              write (nout,"(1x,a,1pe12.4)")    "dx_CartTopog    : ",Partz(Izona)%dx_CartTopog
              write (nout,"(1x,a,1pe12.4)")    "H_res           : ",Partz(Izona)%H_res  
              write (nout,"(1x,a,i12)")        "ID_first_vertex : ",Partz(Izona)%ID_first_vertex
              write (nout,"(1x,a,i12)")        "ID_last_vertex  : ",Partz(Izona)%ID_last_vertex
              write (nout,"(1x,a,i12)")        "plan_reservoir_points: ",Partz(Izona)%plan_reservoir_points
              write (nout,"(1x,a,i12)")        "nag_aux         : ",Partz(Izona)%nag_aux
              do i_point=1,plan_reservoir_points
              write (nout,"(1x,a,3(1pe12.4))") "plan_reservoir_pos   : ",Partz(Izona)%plan_reservoir_pos(i_point,:)                  
              end do
              write (nout,"(1x,a,i12)")        "dam_zone_ID          : ",Partz(Izona)%dam_zone_ID
              write (nout,"(1x,a,i12)")        "dam_zone_n_vertices  : ",Partz(Izona)%dam_zone_n_vertices  
              if (dam_zone_ID>1) then
              do i_point=1,dam_zone_n_vertices
              write (nout,"(1x,a,3(1pe12.4))") "dam_zone_vertices    : ",Partz(Izona)%dam_zone_vertices(i_point,:)                  
              end do
              endif
              endif
!AA504 end

!
            case ("pool")
              Tratto(index)%ColorCode = icolor
              if ( nout > 0 .AND. index == indexi ) write (nout,"(1x,a,z8)")      "Color           : ",Tratto(index)%colorCode

          end select
   
       end do  MULTI_INDEX_LOOP

    end if
!
!.. end of boundary loading
!
    call ReadRiga ( ainp,comment,nrighe,ioerr,ninp )
!
  end do

  return
  end subroutine ReadInputBoundaries
!---split

