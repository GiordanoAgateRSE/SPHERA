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
integer(4) :: n,numv,zone_ID,ipointer,Medium,icolor,icord    
integer(4) :: ioerr,npointv,IC_source_type,Car_top_zone
#ifdef SPACE_3D
integer(4) :: plan_reservoir_points,i_point,dam_zone_ID,dam_zone_n_vertices,i2
integer(4) :: ID_first_vertex_sel,ID_last_vertex_sel
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
   read(token,*,iostat=ioerr) zone_ID
   NumberEntities(8) = max(zone_ID,NumberEntities(8))
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
      ID_first_vertex_sel = 0
      ID_last_vertex_sel = 0
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
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"RRGGBB COLOR",ninp,ulog))   &
            return
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
#ifdef SPACE_3D
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) ID_first_vertex_sel,          &
            ID_last_vertex_sel
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "ID_first_vertex_sel,ID_last_vertex_sel",ninp,ulog)) return
#endif
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
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"FIXED: RRGGBB COLOR",ninp,  &
            ulog)) return
         move = "std"
! Boundary condition "velo"
      case("velo")
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) velocity
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"VELO: NORMAL VELOCITY",ninp,&
            ulog)) return
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
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"VELO: FLOW RATE",ninp,ulog))&
            return
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
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"dx_CartTopog,H_res",  &
                  ninp,ulog)) return
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               if (ioerr==0) read(ainp,*,iostat=ioerr)                         &
                  ID_first_vertex_sel,ID_last_vertex_sel
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "ID_first_vertex_sel,ID_last_vertex_sel",ninp,ulog)) return
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               if (ioerr==0) read(ainp,*,iostat=ioerr) plan_reservoir_points
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "plan_reservoir_points",ninp,ulog)) return
               do i2=1,plan_reservoir_points
                  call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
                  if (ioerr==0) read(ainp,*,iostat=ioerr)                      &
                     plan_reservoir_pos(i2,1),plan_reservoir_pos(i2,2)
                  if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "plan_reservoir_vertices",ninp,ulog)) return
               enddo
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               if (ioerr==0) read(ainp,*,iostat=ioerr) dam_zone_ID,            &
                  dam_zone_n_vertices
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "dam_zone_ID and dam_zone_vertices",ninp,ulog)) return
               if (dam_zone_ID>0) then
                  do i2=1,dam_zone_n_vertices
                     call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
                     if (ioerr==0) read(ainp,*,iostat=ioerr)                   &
                        dam_zone_vertices(i2,1),dam_zone_vertices(i2,2)
                     if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                 &
                        "dam zone vertices",ninp,ulog)) return
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
! Boundary condition "zmax"
      case("zmax")
         NumberEntities(3) = NumberEntities(3) + 1
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Medium
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM INDEX",ninp,ulog))   &
            return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         pressu = trim(GetToken(ainp,1,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"ZMAX PRESSURE TYPE",ninp,   &
            ulog)) return
         if (pressu=="pa") then
! uniform pressure ("pa") is assigned
            token = lcase(GetToken(ainp,2,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,   &
               ulog)) return
            read(token,*,iostat=ioerr) valp
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"PRESSURE VALUES",ninp,   &
               ulog)) return
               elseif (pressu=="qp") then
! hydrostatic pressure based on a reference free surface height ("qp")
                  token = lcase(GetToken(ainp,2,ioerr))
                  if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "INITIAL PIEZO LINE",ninp,ulog)) return
                  read(token,*,iostat=ioerr) valp
                  if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "INITIAL PIEZO LINE",ninp,ulog)) return
                  else
                     if (ulog>0) write(ulog,*) "Unknown option: ",trim(ainp)
                     stop
         endif
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) Car_top_zone
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"Car_top_zone",ninp,ulog))   &
            return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) dx_CartTopog,H_res
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"dx_CartTopog,z_max",ninp,   &
            ulog)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) ID_first_vertex_sel,          &
            ID_last_vertex_sel
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "ID_first_vertex_sel,ID_last_vertex_sel",ninp,ulog)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         if (ioerr==0) read(ainp,*,iostat=ioerr) plan_reservoir_points
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"plan_reservoir_points",ninp,&
            ulog)) return
         do i2=1,plan_reservoir_points
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            if (ioerr==0) read(ainp,*,iostat=ioerr) plan_reservoir_pos(i2,1),  &
               plan_reservoir_pos(i2,2)
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"plan_reservoir_vertices",&
               ninp,ulog)) return
         enddo
#endif
      case default
         write(ulog,*) "Unrecognised BOUNDARY type: ",tipo
         ier = 101
         return
   endselect
   if (input_second_read.eqv..true.) then
      Partz(zone_ID)%label = label
      Partz(zone_ID)%tipo = tipo
      Partz(zone_ID)%Medium = Medium
      Partz(zone_ID)%IC_source_type = IC_source_type
      Partz(zone_ID)%Car_top_zone = Car_top_zone
      Partz(zone_ID)%DBSPH_fictitious_reservoir_flag =                         &
         DBSPH_fictitious_reservoir_flag
#ifdef SPACE_3D
         Partz(zone_ID)%ID_first_vertex_sel = ID_first_vertex_sel
         Partz(zone_ID)%ID_last_vertex_sel = ID_last_vertex_sel
         if ((IC_source_type==2).or.(tipo=="zmax")) then
            Partz(zone_ID)%dx_CartTopog = dx_CartTopog
            Partz(zone_ID)%H_res = H_res
            Partz(zone_ID)%plan_reservoir_points = plan_reservoir_points
            Partz(zone_ID)%plan_reservoir_pos = plan_reservoir_pos
         endif
         if (IC_source_type==2) then
            Partz(zone_ID)%dam_zone_ID = dam_zone_ID
            Partz(zone_ID)%dam_zone_n_vertices = dam_zone_n_vertices
            if (dam_zone_ID>0) Partz(zone_ID)%dam_zone_vertices =              &
               dam_zone_vertices
         endif
#endif
      Partz(zone_ID)%icol = icolor
      Partz(zone_ID)%bend = bends
      Partz(zone_ID)%move = move
      if (npointv/=0) then
         Partz(zone_ID)%npointv = npointv
         Partz(zone_ID)%vlaw(icoordp(0:ncord,ncord-1),1:npointv) =             &
            valuev(0:ncord,1:npointv)
      endif
      Partz(zone_ID)%vel = zero
      Partz(zone_ID)%vel(icoordp(1:ncord,ncord-1)) = values3(1:ncord)
      Partz(zone_ID)%pressure = pressu
      Partz(zone_ID)%valp = valp
      Partz(zone_ID)%slip_coefficient_mode = slip_coefficient_mode
      Partz(zone_ID)%BC_shear_stress_input = BC_shear_stress_input
      if ((pool_plane=="X").or.(pool_plane=="x")) Partz(zone_ID)%ipool = 1
      if ((pool_plane=="Y").or.(pool_plane=="y")) Partz(zone_ID)%ipool = 2
      if ((pool_plane=="Z").or.(pool_plane=="z")) Partz(zone_ID)%ipool = 3
      Partz(zone_ID)%pool = pool_value
      Tratto(zone_ID)%tipo = tipo
#ifdef SPACE_3D
         Tratto(zone_ID)%numvertices = numv
         Tratto(zone_ID)%inivertex = ipointer
#endif
      Tratto(zone_ID)%laminar_no_slip_check = laminar_no_slip_check
      Tratto(zone_ID)%Medium = Medium
      Tratto(zone_ID)%velocity = values1
      Tratto(zone_ID)%NormVelocity = velocity
      Tratto(zone_ID)%FlowRate = flowrate
      Tratto(zone_ID)%zone = zone_ID
      Tratto(zone_ID)%ColorCode = icolor
      if (ulog>0) then
         if (zone_ID>1) write(ulog,*)
         write(ulog,"(1x,a,i5,1x)") "Boundary        : ",zone_ID
         write(ulog,"(1x,a,2x,a)") "Type            : ",Tratto(zone_ID)%tipo 
         if (tipo=="fixe") then
            write(ulog,"(1x,a,l12)") "Laminar no-slip check: ",                &
               Tratto(zone_ID)%laminar_no_slip_check
            elseif (tipo=="peri") then
               write(ulog,"(1x,a,i3)") "Medium Index    : ",                   &
                  Tratto(zone_ID)%Medium
               elseif (tipo=="pool") then
                  write(ulog,"(1x,a,i3)") "Medium Index    : ",                &
                     Tratto(zone_ID)%Medium
                  elseif (tipo=="tapi") then
                     do n=1,ncord
                        icord = icoordp(n,ncord-1)
                        write(ulog,"(1x,a,a,1pe12.4)") xyzlabel(icord),        &
                           " Velocity      : ",Tratto(zone_ID)%velocity(n)
                     enddo
         endif
#ifdef SPACE_2D
            write(ulog,"(1x,a)") "Vertices List"
            write(ulog,"(1x,10i5)")                                            &
BoundaryVertex(Tratto(zone_ID)%inivertex:Tratto(zone_ID)%inivertex+Tratto(zone_ID)%numvertices-1)
            write(ulog,"(1x,a,z6)") "Color           : ",                      &
               Tratto(zone_ID)%colorCode
#endif
      endif
      select case (tipo)
         case("fixe")
#ifdef SPACE_3D
               Tratto(zone_ID)%ColorCode = icolor
               if (ulog>0) write(ulog,"(1x,a,z8)") "Color           : ",       &
                  Tratto(zone_ID)%colorCode
               write(ulog,"(1x,a,i12)") "ID_first_vertex_sel : ",              &
                  Partz(zone_ID)%ID_first_vertex_sel
               write(ulog,"(1x,a,i12)") "ID_last_vertex_sel  : ",              &
                  Partz(zone_ID)%ID_last_vertex_sel
#endif
               write(ulog,"(1x,a,i3)")   "Slip coeff. mode : ",                &
                  Partz(zone_ID)%slip_coefficient_mode
               if (slip_coefficient_mode==1) then
                  write(ulog,"(1x,a,1pe12.4)") "Slip coeff.    : ",            &
                     Partz(zone_ID)%BC_shear_stress_input
                  elseif (slip_coefficient_mode==2) then
                     write(ulog,"(1x,a,1pe12.4)") "Wall mean roughness:",      &
                        Partz(zone_ID)%BC_shear_stress_input
               endif
         case("tapi")
#ifdef SPACE_3D
               Tratto(zone_ID)%ColorCode = icolor
               if (ulog>0) write(ulog,"(1x,a,z8)")                             &
                  "Color           : ",Tratto(zone_ID)%colorCode
#elif defined SPACE_2D
                  numv = Tratto(zone_ID)%numvertices
                  if (numv/=2) then 
                     if (ulog>0) then
                        write(ulog,'(2(a),i15)') "TAPIS boundary type: 2 ",    &
                           "vertices are requested:",numv
                        write(ulog,"(1x,a,1pe12.4)") "Slip coeff.    : ",      &
                           Partz(zone_ID)%BC_shear_stress_input
                     endif
                     ier = 103
                     return
                  endif
#endif
         case("peri")
#ifdef SPACE_2D
               i1 = BoundaryVertex(Tratto(zone_ID)%inivertex)
               i2 =                                                            &
BoundaryVertex(Tratto(zone_ID)%inivertex+Tratto(zone_ID)%numvertices-1)
! Error if the first and the last vertices are different 
               if (i2/=i1) then 
                  if (ulog>0) write(ulog,'(a,2i15)')                           &
"PERIMETER boundary type: first and last vertices are different: ",i2,i1
                  ier = 102
                  return
               endif
               if (ulog>0) then
                  write(ulog,"(1x,a,i3)") "Zone            : ",                &
                     zone_ID,Partz(zone_ID)%label
                  write(ulog,"(1x,a,i3)") "Medium Index    : ",                &
                     Partz(zone_ID)%Medium
                  write(ulog,"(1x,a,Z6.6)") "Color           : ",              &
                     Partz(zone_ID)%icol
                  write(ulog,"(1x,a,2x,a)") "Bends           : ",              &
                     Partz(zone_ID)%bend
                  write(ulog,"(1x,a,2x,a)") "Movement Type   : ",              &
                     Partz(zone_ID)%move
                  if (Partz(zone_ID)%move=="law") then
                     write(ulog,"(1x,2(a),i3)") "Velocity Table - Number of ", &
                        "Points: ",Partz(zone_ID)%npointv
                     do i2=1,Partz(zone_ID)%npointv
                        write(ulog,"(a,i3,1p,4(2x,a,e12.4))") " Point",i2,     &
                           (xyzlabel(icoordp(n,ncord-1)),                      &
                           Partz(zone_ID)%vlaw(icoordp(n,ncord-1),i2),n=0,ncord)
                     enddo
                  endif
                  do n=1,ncord
                     icord = icoordp(n,ncord-1)
                     write(ulog,"(1x,a,a,1pe12.4)") xyzlabel(icord),           &
                        " velocity       : ",Partz(zone_ID)%vel(icord) 
                  enddo
                  write(ulog,"(1x,a,2x,a)")    "Pressure Type   : ",           &
                     Partz(zone_ID)%pressure
                  write(ulog,"(1x,a,1pe12.4)") "Pressure Value  : ",           &
                     Partz(zone_ID)%valp
               endif       
#elif defined SPACE_3D
                  Tratto(zone_ID)%ColorCode = icolor
                  if (ulog>0) write(ulog,"(1x,a,z8)") "Color           : ",    &
                     Tratto(zone_ID)%colorCode
#endif
            write(ulog,"(1x,a,i12)") "IC_source_type  : ",                     &
               Partz(zone_ID)%IC_source_type
            write(ulog,"(1x,a,i12)") "Car_top_zone    : ",                     &
               Partz(zone_ID)%Car_top_zone
            write(ulog,"(1x,a,l12)") "DBSPH_fictitious_reservoir_flag : ",     &
               Partz(zone_ID)%DBSPH_fictitious_reservoir_flag
#ifdef SPACE_3D
               if (IC_source_type==2) then
                  write(ulog,"(1x,a,1pe12.4)") "dx_CartTopog    : ",           &
                     Partz(zone_ID)%dx_CartTopog
                  write(ulog,"(1x,a,1pe12.4)") "H_res           : ",           &
                     Partz(zone_ID)%H_res  
                  write(ulog,"(1x,a,i12)") "ID_first_vertex_sel : ",           &
                     Partz(zone_ID)%ID_first_vertex_sel
                  write(ulog,"(1x,a,i12)") "ID_last_vertex_sel  : ",           &
                     Partz(zone_ID)%ID_last_vertex_sel
                  write(ulog,"(1x,a,i12)") "plan_reservoir_points: ",          &
                     Partz(zone_ID)%plan_reservoir_points
                  do i_point=1,plan_reservoir_points
                     write(ulog,"(1x,a,3(1pe12.4))") "plan_reservoir_pos   : ",&
                        Partz(zone_ID)%plan_reservoir_pos(i_point,:)
                  enddo
                  write(ulog,"(1x,a,i12)") "dam_zone_ID          : ",          &
                     Partz(zone_ID)%dam_zone_ID
                  write(ulog,"(1x,a,i12)") "dam_zone_n_vertices  : ",          &
                     Partz(zone_ID)%dam_zone_n_vertices  
                  if (dam_zone_ID>0) then
                     do i_point=1,dam_zone_n_vertices
                        write(ulog,"(1x,2(a),3(1pe12.4))") "dam_zone_vertices",&
                           "    : ",Partz(zone_ID)%dam_zone_vertices(i_point,:)                  
                     enddo
                  endif
               endif
         case("zmax")
            write(ulog,"(1x,a,i3)") "Medium Index    : ",                      &
               Partz(zone_ID)%Medium
            write(ulog,"(1x,a,2x,a)") "Pressure Type   : ",                    &
               Partz(zone_ID)%pressure
            write(ulog,"(1x,a,1pe12.4)") "Pressure Value  : ",                 &
               Partz(zone_ID)%valp
            write(ulog,"(1x,a,i12)") "Car_top_zone    : ",                     &
               Partz(zone_ID)%Car_top_zone
            write(ulog,"(1x,a,1pe12.4)") "dx_CartTopog    : ",                 &
               Partz(zone_ID)%dx_CartTopog
            write(ulog,"(1x,a,1pe12.4)") "z_max           : ",                 &
               Partz(zone_ID)%H_res
            write(ulog,"(1x,a,i12)") "ID_first_vertex_sel : ",                 &
               Partz(zone_ID)%ID_first_vertex_sel
            write(ulog,"(1x,a,i12)") "ID_last_vertex_sel  : ",                 &
               Partz(zone_ID)%ID_last_vertex_sel
            write(ulog,"(1x,a,i12)") "plan_zmax_zone_points: ",                &
               Partz(zone_ID)%plan_reservoir_points
            do i_point=1,plan_reservoir_points
               write(ulog,"(1x,a,3(1pe12.4))") "plan_zmax_zone_pos   : ",      &
                  Partz(zone_ID)%plan_reservoir_pos(i_point,:)
            enddo
#endif
         case("pool")
            Tratto(zone_ID)%ColorCode = icolor
            if (ulog>0) write(ulog,"(1x,a,z8)") "Color           : ",          &
               Tratto(zone_ID)%colorCode
      endselect
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputBoundaries
