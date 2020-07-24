!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2020 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.9.0.0
! SPHERA v.9.0.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: ReadInputBoundaries                    
! Description: Reading input data for the boundary treatment scheme SA-SPH 
!              (semi-analytic approach; Di Monaco et al., 2011, EACFM).                      
!-------------------------------------------------------------------------------
subroutine ReadInputBoundaries(NumberEntities,Partz,Tratto,                    &
#ifdef SPACE_2D
   BoundaryVertex,                                                             &
#endif
   ainp,comment,nrighe,ier,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: DBSPH_fictitious_reservoir_flag,laminar_no_slip_check
integer(4) :: nrighe,ier,ninp,ulog,slip_coefficient_mode
integer(4),dimension(20) :: NumberEntities
#ifdef SPACE_2D
integer(4),dimension(NumBVertices) :: BoundaryVertex
#endif
character(1) :: comment
character(len=lencard) :: ainp
type (TyZone),dimension(NPartZone) :: Partz
type (TyBoundaryStretch),dimension(NumTratti) :: Tratto
integer(4) :: n,index,numv,indexi,indexf,Izona,ipointer,Medium,icolor,icord    
integer(4) :: ioerr,npointv,IC_source_type,Car_top_zone
#ifdef SPACE_3D
integer(4) :: plan_reservoir_points,i_point,dam_zone_ID,dam_zone_n_vertices,i2
integer(4) :: ID_first_vertex,ID_last_vertex
#elif defined SPACE_2D
integer(4) :: i2,i1
#endif
double precision :: pool_value,velocity,valp,flowrate,BC_shear_stress_input
#ifdef SPACE_3D
double precision :: dx_CartTopog,H_res
#endif
double precision,dimension(3) :: values1,values3
double precision,dimension(0:3,maxpointsvlaw) :: valuev
#ifdef SPACE_3D
double precision :: plan_reservoir_pos(4,2),dam_zone_vertices(4,2)
#endif
character(1) :: pool_plane,bends,slip
character(2) :: pressu
character(3) :: move
character(4) :: tipo
character(6) :: token_color
character(8) :: label
character(100) :: token
logical,external :: ReadCheck
integer(4),external :: ptcolorrgb
character(100), external :: lcase,GetToken
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
npointv = 0
values3 = zero
valp = zero
slip_coefficient_mode = 0
BC_shear_stress_input = -9.99d9
!------------------------
! Statements
!------------------------
! In case of restart, input data are not read
if (restart) then
   do while (trim(lcase(ainp))/="##### end boundaries #####")
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BOUNDARIES DATA",ninp,ulog))   &
         return
   enddo
   return
endif
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BOUNDARIES DATA",ninp,ulog)) return
! Reading input data
do while (trim(lcase(ainp))/="##### end boundaries #####")
   label = ainp(1:8)
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BOUNDARIES INDEX",ninp,ulog))     &
      return
   token = GetToken(ainp,1,ioerr)
   read(token,*,iostat=ioerr) indexi
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BOUNDARY INDEX",ninp,ulog)) return
   indexf = indexi
   token = GetToken(ainp,2,ioerr)
   if (token/="") read(token,*,iostat=ioerr) indexf
   NumberEntities(8) = max(indexf,NumberEntities(8))
! Boundary type
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"BOUNDARY TYPE",ninp,ulog)) return
   tipo = trim(lcase(ainp(1:4)))
   numv = 0
   ipointer = 0
   move = "   "
   Medium = 0
   icolor = Z'FF0000'  
   values1 = zero
   values3 = zero
   pool_plane = " "
   pool_value = zero
   slip_coefficient_mode = 0
   velocity = zero
   flowrate = zero
   pressu = "  "
   valp = zero
   IC_source_type = 0
   Car_top_zone = 0
   DBSPH_fictitious_reservoir_flag = .false.
#ifdef SPACE_3D
      dx_CartTopog = 0.d0
      H_res = 0.d0
      ID_first_vertex = 0
      ID_last_vertex = 0
      plan_reservoir_points = 0
      dam_zone_ID = 0
      dam_zone_n_vertices = 0
      plan_reservoir_pos = 0.d0 
      dam_zone_vertices = 0.d0
#endif
   select case (tipo)
! Boundary condition "leve", "crit" or "open"
      case("leve","crit","open")    
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         token = GetToken(ainp,1,ioerr)
         token_color(1:2) = token(5:6)
         token_color(3:4) = token(3:4)
         token_color(5:6) = token(1:2) 
         read(token_color,'(Z6)',iostat=ioerr) icolor
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,  &
            ulog)) return
! Boundary condition "fixe"
      case("fixe")    
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) slip_coefficient_mode,        &
            laminar_no_slip_check
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "FIXED: SLIP COEFFICIENT MODE, LAMINAR NO-SLIP CHECK",ninp,ulog))  &
            return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) BC_shear_stress_input
         if (slip_coefficient_mode==1) then
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SLIP COEFFICIENT",ninp,  &
               ulog)) return
            elseif (slip_coefficient_mode==2) then
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"WALL MEAN ROUGHNESS", &
                  ninp,ulog)) return
         endif
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         token = GetToken(ainp,1,ioerr)
         token_color(1:2) = token(5:6)
         token_color(3:4) = token(3:4)
         token_color(5:6) = token(1:2) 
         read(token_color,'(Z6)',iostat=ioerr) icolor
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,  &
            ulog)) return
! Boundary condition "sour"
      case("sour")    
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) Medium
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: MEDIUM INDEX",ninp, &
            ulog)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE ",ninp,   &
            ulog)) return
         token = GetToken(ainp,1,ioerr)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE",ninp,    &
            ulog)) return
         read(token,*,iostat=ioerr) flowrate
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE: FLOW RATE",ninp,    &
            ulog)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         pressu = trim(GetToken(ainp,1,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE TYPE",ninp, &
            ulog)) return
         if (pressu=="pa") then  
            token = GetToken(ainp,2,ioerr)
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE VALUES", &
               ninp,ulog)) return
            read(token,*,iostat=ioerr) valp
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PRESSURE VALUES", &
               ninp,ulog)) return
            elseif (pressu=="qp") then
               token = GetToken(ainp,2,ioerr)
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PIEZO LINE",   &
                  ninp,ulog)) return
               read(token,*,iostat=ioerr) valp
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"SOURCE PIEZO LINE",   &
                  ninp,ulog)) return
               else
                  if (ulog>0) write(ulog,*) "Unknown option: ",trim(ainp),     &
                     " in source boundary."
                  stop
         endif
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         token = GetToken(ainp,1,ioerr)
         token_color(1:2) = token(5:6)
         token_color(3:4) = token(3:4)
         token_color(5:6) = token(1:2) 
         read(token_color,'(Z6)',iostat=ioerr) icolor
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp&
            ,ulog)) return
         move   = "std"
! Boundary condition "velo"
      case("velo")
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) velocity
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "VELO: NORMAL VELOCITY",ninp,ulog)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         token = GetToken(ainp,1,ioerr)
         token_color(1:2) = token(5:6)
         token_color(3:4) = token(3:4)
         token_color(5:6) = token(1:2) 
         read(token_color,'(Z6)',iostat=ioerr) icolor
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,  &
            ulog)) return
         move = "std"
         pressu = "pa"
         valp = zero
! Boundary condition "flow"
      case("flow")    
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) flowrate
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"VELO: FLOW RATE",ninp,ulog) &
            ) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         token = GetToken(ainp,1,ioerr)
         token_color(1:2) = token(5:6)
         token_color(3:4) = token(3:4)
         token_color(5:6) = token(1:2) 
         read(token_color,'(Z6)',iostat=ioerr) icolor
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,  &
            ulog)) return
         move = "std"
         pressu = "pa"
         valp = zero
! Boundary condition "peri"
      case("peri")    
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadInputParticlesData(NumberEntities,Medium,icolor,bends,move,  &
            slip,npointv,valuev,values3,pressu,valp,ainp,comment,nrighe,ier,   &
            ninp,ulog)
         if (ier/=0) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) IC_source_type,Car_top_zone,  &
            DBSPH_fictitious_reservoir_flag
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "IC_source_type, Car_top_zone, DBSPH_fictitious_reservoir_flag",   &
            ninp,ulog)) return
#ifdef SPACE_3D
         if (IC_source_type==2) then
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            if (ioerr==0) read(ainp,*,iostat=ioerr) dx_CartTopog,H_res
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"dx_CartTopog,H_res",ninp,&
               ulog)) return
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            if (ioerr==0) read(ainp,*,iostat=ioerr)                            &
               ID_first_vertex,ID_last_vertex
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "ID_first_vertex,ID_last_vertex",ninp,ulog)) return
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            if (ioerr==0) read(ainp,*,iostat=ioerr) plan_reservoir_points
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "plan_reservoir_points",ninp,ulog)) return
            do i2=1,plan_reservoir_points
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               if (ioerr==0) read(ainp,*,iostat=ioerr)                         &
                  plan_reservoir_pos(i2,1),plan_reservoir_pos(i2,2)
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "plan_reservoir_vertices",ninp,ulog)) return
            enddo
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            if (ioerr==0) read(ainp,*,iostat=ioerr) dam_zone_ID,               &
               dam_zone_n_vertices
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "dam_zone_ID and dam_zone_vertices",ninp,ulog)) return
            if (dam_zone_ID>0) then
               do i2=1,dam_zone_n_vertices
                  call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
                  if (ioerr==0) read(ainp,*,iostat=ioerr)                      &
                     dam_zone_vertices(i2,1),dam_zone_vertices(i2,2)
                  if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"dam zone vertices",&
                     ninp,ulog)) return
               enddo
            endif
         endif
#endif
! Boundary condition "tapi"
      case("tapi")
#ifdef SPACE_3D
! It returns an error if the number of vertices is not 2
            if (numv/=2) then 
               write(ulog,'(a,i15)')                                           &
                  "TAPIS boundary type: 2 vertices are requested:",numv
               ier = 103
               return
            endif
#endif
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) BC_shear_stress_input
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "TAPIS: SLIP COEFFICIENT",ninp,ulog)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: VELOCITY COMPONENTS",&
            ninp,ulog)) return
         do n=1,ncord
            icord = icoordp(n,ncord-1)
            token = GetToken(ainp,n,ioerr)
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               xyzlabel(icord)//" VELOCITY COMPONENT (TAPIS)",ninp,ulog)) return
            read(token,*,iostat=ioerr) values1(icord)
         enddo
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         token = GetToken(ainp,1,ioerr)
         token_color(1:2) = token(5:6)
         token_color(3:4) = token(3:4)
         token_color(5:6) = token(1:2) 
         read(token_color,'(Z6)',iostat=ioerr) icolor
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"TAPIS: RRGGBB COLOR",ninp,  &
            ulog)) return
! Boundary condition "pool" (only in 3D)
      case("pool")
#ifdef SPACE_3D
            NumberEntities(3) = NumberEntities(3) + 1
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            pool_plane = trim(GetToken(ainp,1,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"POOL: X/Y/Z/ PLANE LABEL"&
               ,ninp,ulog)) return
            token = GetToken(ainp,2,ioerr)
            if (ioerr==0) read(token,*,iostat=ioerr) pool_value
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"POOL: PLANE VALUE",ninp, &
               ulog)) return
            call ReadInputParticlesData(NumberEntities,Medium,icolor,bends,    &
                                        move,slip,npointv,valuev,values3,      &
                                        pressu,valp,ainp,comment,nrighe,ier,   &
                                        ninp,ulog)
            if (ier/=0) return  
#endif 
      case default
         write(ulog,*) "Unrecognised BOUNDARY type: ",tipo
         ier = 101
         return
   end  select
   if (input_second_read.eqv..true.) then
      Izona = NumberEntities(3)
      Partz(Izona)%label = label
      Partz(Izona)%tipo = tipo
      Partz(Izona)%Medium = Medium
      Partz(Izona)%IC_source_type = IC_source_type
      Partz(Izona)%Car_top_zone = Car_top_zone
      Partz(Izona)%DBSPH_fictitious_reservoir_flag =                           &
         DBSPH_fictitious_reservoir_flag
#ifdef SPACE_3D
      if (IC_source_type==2) then
         Partz(Izona)%dx_CartTopog = dx_CartTopog
         Partz(Izona)%H_res = H_res
         Partz(Izona)%ID_first_vertex = ID_first_vertex
         Partz(Izona)%ID_last_vertex = ID_last_vertex
         Partz(Izona)%plan_reservoir_points = plan_reservoir_points
         Partz(Izona)%plan_reservoir_pos = plan_reservoir_pos
         Partz(Izona)%dam_zone_ID = dam_zone_ID
         Partz(Izona)%dam_zone_n_vertices = dam_zone_n_vertices
         if (dam_zone_ID>0) Partz(Izona)%dam_zone_vertices = dam_zone_vertices
      endif
#endif
      Partz(Izona)%icol = icolor
      Partz(Izona)%bend = bends
      Partz(Izona)%move = move
      if (npointv/=0) then
         Partz(Izona)%npointv = npointv
         Partz(Izona)%vlaw(icoordp(0:ncord,ncord-1),1:npointv) =               &
            valuev(0:ncord,1:npointv)
      endif
      Partz(Izona)%vel = zero
      Partz(Izona)%vel(icoordp(1:ncord,ncord-1)) = values3(1:ncord)
      Partz(Izona)%pressure = pressu
      Partz(Izona)%valp = valp
      Partz(Izona)%Indix(1) = indexi
      Partz(Izona)%Indix(2) = indexf
      Partz(Izona)%slip_coefficient_mode = slip_coefficient_mode
      Partz(Izona)%BC_shear_stress_input = BC_shear_stress_input
      if ((pool_plane=="X").or.(pool_plane=="x")) Partz(Izona)%ipool = 1
      if ((pool_plane=="Y").or.(pool_plane=="y")) Partz(Izona)%ipool = 2
      if ((pool_plane=="Z").or.(pool_plane=="z")) Partz(Izona)%ipool = 3
      Partz(Izona)%pool = pool_value
! Constraints
      MULTI_INDEX_LOOP: do index=indexi,indexf       
         Tratto(index)%tipo = tipo
#ifdef SPACE_3D
            Tratto(index)%numvertices = numv
            Tratto(index)%inivertex = ipointer
#endif
         Tratto(index)%laminar_no_slip_check = laminar_no_slip_check
         Tratto(index)%Medium = Medium
         Tratto(index)%velocity = values1
         Tratto(index)%NormVelocity = velocity
         Tratto(index)%FlowRate = flowrate
         Tratto(index)%zone = Izona
         Tratto(index)%ColorCode = icolor
         if ((ulog>0).and.(index==indexi)) then
            if (index>1) write(ulog,*)
            if (indexf==indexi) write(ulog,"(1x,a,i5,1x,a)")                   &
               "Boundary        : ",indexi
            if (indexf/=indexi) write(ulog,"(1x,a,i5,1x,a,i5)")                &
               "Boundary        : ",indexi,"   to",indexf
            write(ulog,"(1x,a,2x,a)") "Type            : ",Tratto(index)%tipo 
            if (tipo=="fixe") then
               write(ulog,"(1x,a,l12)") "Laminar no-slip check: ",             &
                  Tratto(index)%laminar_no_slip_check
               elseif (tipo=="peri") then
                  write(ulog,"(1x,a,i3)") "Medium Index    : ",                &
                     Tratto(index)%Medium
                  elseif (tipo=="pool") then
                     write(ulog,"(1x,a,i3,1x,a)") "Medium Index    : ",        &
                        Tratto(index)%Medium
                     elseif (tipo=="tapi") then
                        do n=1,ncord
                           icord = icoordp(n,ncord-1)
                           write(ulog,"(1x,a,a,1pe12.4)") xyzlabel(icord),     &
                              " Velocity      : ",Tratto(index)%velocity(n)
                        enddo
            endif
#ifdef SPACE_2D
               write(ulog,"(1x,a)") "Vertices List"
               write(ulog,"(1x,10i5)")                                         &
BoundaryVertex(Tratto(index)%inivertex:Tratto(index)%inivertex+Tratto(index)%numvertices-1)
               write(ulog,"(1x,a,z6)") "Color           : ",Tratto(index)%colorCode
#endif
         endif
         select case (tipo)
            case("fixe")
#ifdef SPACE_3D
                  Tratto(index)%ColorCode = icolor
                  if ((ulog>0).and.(index==indexi)) write(ulog,                &
                     "(1x,a,z8)") "Color           : ",                        &
                     Tratto(index)%colorCode
                  Tratto(index)%ColorCode = icolor
#endif
                  write(ulog,"(1x,a,i3)")   "Slip coeff. mode : ",             &
                     Partz(Izona)%slip_coefficient_mode
                  if (slip_coefficient_mode==1) then
                     write(ulog,"(1x,a,1pe12.4)") "Slip coeff.    : ",         &
                        Partz(Izona)%BC_shear_stress_input
                     elseif (slip_coefficient_mode==2) then
                        write(ulog,"(1x,a,1pe12.4)") "Wall mean roughness:",   &
                           Partz(Izona)%BC_shear_stress_input
                  endif
            case("tapi")
#ifdef SPACE_3D
                  Tratto(index)%ColorCode = icolor
                  if (ulog>0.and.index==indexi) write(ulog,"(1x,a,z8)")        &
                     "Color           : ",Tratto(index)%colorCode
                  Tratto(index)%ColorCode = icolor
#elif defined SPACE_2D
                     numv = Tratto(index)%numvertices
                     if (numv/=2) then 
                        if (ulog>0) then
                           write(ulog,'(a,i15)')                               &
                              "TAPIS boundary type: 2 vertices are requested:",&
                              numv
                           write(ulog,"(1x,a,1pe12.4)") "Slip coeff.    : ",   &
                              Partz(Izona)%BC_shear_stress_input
                        endif
                        ier = 103
                        return
                     endif
#endif
            case("peri")
#ifdef SPACE_2D
                  i1 = BoundaryVertex(Tratto(index)%inivertex)
                  i2 =                                                         &
BoundaryVertex(Tratto(index)%inivertex+Tratto(index)%numvertices-1)
! Error if the first and the last vertices are different 
                  if (i2/=i1) then 
                     if (ulog>0) write(ulog,'(a,2i15)')                        &
"PERIMETER boundary type: first and last vertices are different: ",i2,i1
                     ier = 102
                     return
                  endif
                  if (ulog>0) then
                     write(ulog,"(1x,a,i3,1x,a)")                              &
                        "Zone            : ",Izona,Partz(Izona)%label
                     write(ulog,"(1x,a,i3)") "Medium Index    : ",             &
                        Partz(Izona)%Medium
                     write(ulog,"(1x,a,Z6.6)") "Color           : ",           &
                        Partz(Izona)%icol
                     write(ulog,"(1x,a,2x,a)") "Bends           : ",           &
                        Partz(Izona)%bend
                     write(ulog,"(1x,a,2x,a)") "Movement Type   : ",           &
                        Partz(Izona)%move
                     if (Partz(Izona)%move=="law") then
                        write(ulog,"(1x,a,i3)")                                & 
                           "Velocity Table - Number of Points: ",              &
                           Partz(Izona)%npointv
                        do i2=1,Partz(Izona)%npointv
                           write(ulog,"(a,i3,1p,4(2x,a,e12.4))") " Point",i2,  &
                              (xyzlabel(icoordp(n,ncord-1)),                   &
                              Partz(Izona)%vlaw(icoordp(n,ncord-1),i2),n=0,    &
                              ncord)
                        enddo
                     endif
                     do n=1,ncord
                        icord = icoordp(n,ncord-1)
                        write(ulog,"(1x,a,a,1pe12.4)") xyzlabel(icord),        &
                           " velocity       : ",Partz(Izona)%vel(icord) 
                     enddo
                     write(ulog,"(1x,a,2x,a)")    "Pressure Type   : ",        &
                        Partz(Izona)%pressure
                     write(ulog,"(1x,a,1pe12.4)") "Pressure Value  : ",        &
                        Partz(Izona)%valp
                  endif       
#elif defined SPACE_3D
                     Tratto(index)%ColorCode = icolor
                     if ((ulog>0).and.(index==indexi)) write(ulog,"(1x,a,z8)") &
                        "Color           : ",Tratto(index)%colorCode
#endif
               write(ulog,"(1x,a,i12)") "IC_source_type  : ",                  &
                  Partz(Izona)%IC_source_type
               write(ulog,"(1x,a,i12)") "Car_top_zone    : ",                  &
                  Partz(Izona)%Car_top_zone
               write(ulog,"(1x,a,l12)") "DBSPH_fictitious_reservoir_flag : ",  &
                  Partz(Izona)%DBSPH_fictitious_reservoir_flag
#ifdef SPACE_3D
               if (IC_source_type==2) then
                  write(ulog,"(1x,a,1pe12.4)") "dx_CartTopog    : ",           &
                     Partz(Izona)%dx_CartTopog
                  write(ulog,"(1x,a,1pe12.4)") "H_res           : ",           &
                     Partz(Izona)%H_res  
                  write(ulog,"(1x,a,i12)") "ID_first_vertex : ",               &
                     Partz(Izona)%ID_first_vertex
                  write(ulog,"(1x,a,i12)") "ID_last_vertex  : ",               &
                     Partz(Izona)%ID_last_vertex
                  write(ulog,"(1x,a,i12)") "plan_reservoir_points: ",          &
                     Partz(Izona)%plan_reservoir_points
                  do i_point=1,plan_reservoir_points
                     write(ulog,"(1x,a,3(1pe12.4))") "plan_reservoir_pos   : " &
                        ,Partz(Izona)%plan_reservoir_pos(i_point,:)                  
                  enddo
                  write(ulog,"(1x,a,i12)") "dam_zone_ID          : ",          &
                     Partz(Izona)%dam_zone_ID
                  write(ulog,"(1x,a,i12)") "dam_zone_n_vertices  : ",          &
                     Partz(Izona)%dam_zone_n_vertices  
                  if (dam_zone_ID>0) then
                     do i_point=1,dam_zone_n_vertices
                        write(ulog,"(1x,a,3(1pe12.4))")                        &
                           "dam_zone_vertices    : ",                          &
                           Partz(Izona)%dam_zone_vertices(i_point,:)                  
                     enddo
                  endif
               endif
#endif
            case("pool")
               Tratto(index)%ColorCode = icolor
               if ((ulog>0).and.(index==indexi)) write(ulog,"(1x,a,z8)")       &
                  "Color           : ",Tratto(index)%colorCode
         endselect
      enddo MULTI_INDEX_LOOP
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputBoundaries
