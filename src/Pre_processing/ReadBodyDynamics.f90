!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: ReadBodyDynamics                    
! Description: Reading input data for body trasnport in fluid flows (Amicarelli et al., 2015, CAF).                  
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ReadBodyDynamics (ainp,comment,nrighe,ier,ninp,nout)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module                            
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,ninp,nout,ioerr,i,Id_body,n_elem,j,Id_elem
integer(4) :: imposed_kinematics,n_records,Ic_imposed 
double precision :: mass
integer(4) :: normal_act(6)
double precision :: L_geom(3),x_CM(3),alfa(3),u_CM(3),omega(3),x_rotC(3)
double precision :: mass_deact(6)
double precision :: Ic(3,3)
character(1) :: comment
character(80) :: ainp,lcase !,token,GetToken
logical,external :: ReadCheck
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
! In case of restart, input data are not read
if (restart) then
   do while (TRIM(lcase(ainp)) /= "##### end body dynamics #####")
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout))&
         return
   enddo
   return
endif
call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout))      &
   return
do while (TRIM(lcase(ainp)) /= "##### end body dynamics #####")
! Reading the number of bodies and the ratio between fluid and body particle 
! size
   read(ainp,*,iostat=ioerr) n_bodies,dx_dxbodies,imping_body_grav
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS GENERAL INPUT",ninp,&
      nout)) return
! Writing the number of bodies and "dx_dxbodies" on the log file
   if ((ncord>0).and.(nout>0)) then
      write(nout,"(1x,a,1p,i12)") "n_bodies:.....................",n_bodies
      write(nout,"(1x,a,1p,e12.4)") "dx_dxbodies:..................",          &
         dx_dxbodies   
      write(nout,"(1x,a,1p,i12)") "imping_body_grav:.............",            &
         imping_body_grav
      write(nout,"(1x,a)")  " "
   endif
! Allocation of the array of the bodies
   if (allocated(body_arr)) then
      else
         allocate(body_arr(n_bodies))  
   endif
! Loop over the transported bodies
   do i=1,n_bodies
! Reading the body parameters
      call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) Id_body,n_elem
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_BODY-N_ELEM",ninp,nout))    &  
         return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) mass
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MASS",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) Ic_imposed
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC_IMPOSED",ninp,nout)) return
      if (Ic_imposed==1) then
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Ic(1,1),Ic(1,2),Ic(1,3)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(1,1-3)",ninp,nout))      &
            return
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Ic(2,1),Ic(2,2),Ic(2,3)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(2,1-3)",ninp,nout))      &
            return
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Ic(3,1),Ic(3,2),Ic(3,3)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(3,1-3)",ninp,nout))      &
            return
      endif
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) alfa(1),alfa(2),alfa(3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALFA",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) x_rotC(1),x_rotC(2),x_rotC(3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_ROTC",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) u_CM(1),u_CM(2),u_CM(3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"U_CM",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) omega(1),omega(2),omega(3)
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"OMEGA",ninp,nout)) return
      call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
      read(ainp,*,iostat=ioerr) imposed_kinematics,n_records
      if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS",ninp,nout))   &
         return       
! Assignment of the body parameters 
      body_arr(Id_body)%n_elem = n_elem
      body_arr(Id_body)%mass = mass
      body_arr(Id_body)%x_CM = x_CM
      body_arr(Id_body)%Ic_imposed = Ic_imposed
      if (Ic_imposed==1) then
         body_arr(Id_body)%Ic = Ic
         else
            body_arr(Id_body)%Ic = 0.d0
      endif   
      body_arr(Id_body)%alfa = alfa
      body_arr(Id_body)%x_rotC = x_rotC
      body_arr(Id_body)%u_CM = u_CM
      body_arr(Id_body)%omega = omega
      body_arr(Id_body)%imposed_kinematics = imposed_kinematics
      body_arr(Id_body)%n_records = n_records
! Writing on the log file
      if (ncord>0) then
         if (nout>0) then
            write(nout,"(1x,a,1p,i12)") "body:.......................",        &
               Id_body
            write(nout,"(1x,a,1p,e12.4)") "mass:.......................",mass
            write(nout,"(1x,a,1p,3e12.4)") "x_CM:.......................",x_CM
            write(nout,"(1x,a,1p,i12)") "IC_imposed:.................",        &
               Ic_imposed
            if (Ic_imposed==1) then
               write(nout,"(1x,a,1p,3e12.4)") "Ic(1,1-3):..................",  &
                  Ic(1,1),Ic(1,2),Ic(1,3)
               write(nout,"(1x,a,1p,3e12.4)") "Ic(2,1-3):..................",  &
                  Ic(2,1),Ic(2,2),Ic(2,3)
               write(nout,"(1x,a,1p,3e12.4)") "Ic(3,1-3):..................",  &
                  Ic(3,1),Ic(3,2),Ic(3,3)
            endif
            write(nout,"(1x,a,1p,3e12.4)") "alfa:.......................",     &
               alfa(1),alfa(2),alfa(3) 
            write(nout,"(1x,a,1p,3e12.4)") "x_rotC:.....................",     &
               x_rotC 
            write(nout,"(1x,a,1p,3e12.4)") "u_CM:.......................",     &
               u_CM
            write(nout,"(1x,a,1p,3e12.4)") "omega:......................",     &
               omega
            write(nout,"(1x,a,1p,i12)") "imposed_kinematics:.........",        &
               imposed_kinematics    
            write(nout,"(1x,a,1p,i12)") "n_records:..................",        &
               n_records
            write(nout,"(1x,a)")  " "
         endif
      endif
! Allocating body elements
      if (ncord>0) then
         else
            allocate(body_arr(Id_body)%elem(body_arr(Id_body)%n_elem))
            allocate(body_arr(Id_body)%body_kinematics(n_records,7))
      endif
! Reading the eventual imposed kinematics
      do j=1,n_records
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) body_arr(Id_body)%body_kinematics(j,:) 
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS_RECORDS",   &
            ninp,nout)) return            
      enddo       
! Reading element parameters                     
      do j=1,n_elem
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Id_elem
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_ELEM",ninp,nout)) return
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) L_geom(1),L_geom(2),L_geom(3)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"L_GEOM",ninp,nout)) return
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM_ELEM",ninp,nout))      &
            return
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) alfa(1),alfa(2),alfa(3)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALFA_ELEM",ninp,nout))      &
            return
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) normal_act(1),normal_act(2),normal_act(3),  &
            normal_act(4),normal_act(5),normal_act(6)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NORMAL_ACT",ninp,nout))     &
            return
         call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) mass_deact(1),mass_deact(2),mass_deact(3),  &
            mass_deact(4),mass_deact(5),mass_deact(6)
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MASS_DEACT",ninp,nout))     &
            return
! Assignment of the element parameters
         body_arr(Id_body)%elem(Id_elem)%L_geom = L_geom
         body_arr(Id_body)%elem(Id_elem)%x_CM = x_CM
         body_arr(Id_body)%elem(Id_elem)%alfa = alfa
         body_arr(Id_body)%elem(Id_elem)%normal_act = normal_act
         body_arr(Id_body)%elem(Id_elem)%mass_deact = mass_deact          
! Writing on the log file
         if (ncord>0) then
            if (nout>0) then
               write(nout,"(1x,a,1p,i12)") "element:....................",     &
                  Id_elem
               write(nout,"(1x,a,1p,3e12.4)") "L_geom:.....................",  &
                  L_geom
               write(nout,"(1x,a,1p,3e12.4)") "x_CM_elem:..................",  &
                  x_CM
               write(nout,"(1x,a,1p,3e12.4)") "alfa_elem:..................",  &
                  alfa(1),alfa(2),alfa(3)    
               write(nout,"(1x,a,1p,6i12)") "normal_act:.................",    &
                  normal_act    
               write(nout,"(1x,a,1p,6e12.4)") "mass_deact:.................",  &
                  mass_deact                                          
               write(nout,"(1x,a)")  " "
            endif
         endif          
      enddo
   enddo         
   call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout))   &
      return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadBodyDynamics

