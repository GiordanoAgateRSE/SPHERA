!AA501b (whole subroutine)
!cfile ReadBodyDynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ReadInputFile.f90
!
! Creation      : 13nov12 Amicarelli-Agate
!
!************************************************************************************
! Module purpose : Reading input for body dynamics
!
! Calling routines: ReadInput
!
! Called routines: /
!
!************************************************************************************
subroutine ReadBodyDynamics (ainp,comment,nrighe,ier,ninp,nout)

! modules
use GLOBAL_MODULE                            
use AdM_USER_TYPE
use ALLOC_Module

!Declarations
implicit none

integer(4)    :: nrighe,ier,ninp,nout,ioerr,i,Id_body,n_elem,j,Id_elem,imposed_kinematics,n_records,Ic_imposed !,k
character( 1) :: comment
character(80) :: ainp,lcase !,token,GetToken
double precision :: mass                          
double precision :: L_geom(3),x_CM(3),alfa(3),u_CM(3),omega(3),mass_deact(6),Ic(3,3),x_rotC(3)
integer(4) :: normal_act(6)

logical,       external :: ReadCheck

! in case of restart the cards are not read
if (restart) then
  do while (TRIM(lcase(ainp)) /= "##### end body dynamics #####")
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout)) return
  end do
  return
end if

call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout)) return

 do while (TRIM(lcase(ainp)) /= "##### end body dynamics #####")
!Reading the number of bodies and the ratio between fluid and body particle size
    read (ainp,*,iostat=ioerr) n_bodies,dx_dxbodies,imping_body_grav
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS GENERAL INPUT",ninp,nout)) return
! Writing the number of bodies and dx_dxbodies on the log file
!AA504 sub
    if ((ncord>0).and.(nout > 0)) then
       write (nout,"(1x,a,1p,i12)")   "n_bodies:.....................",n_bodies
       write (nout,"(1x,a,1p,e12.4)") "dx_dxbodies:..................",dx_dxbodies   
       write (nout,"(1x,a,1p,i12)")   "imping_body_grav:.............",imping_body_grav
       write (nout,"(1x,a)")  " "
    end if
! Allocation of the array of the bodies
    if (allocated(body_arr)) then
    else
    allocate(body_arr(n_bodies))  
    endif
! Loop over the transported bodies
    do i=1,n_bodies
!Reading the body parameters
       call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) Id_body,n_elem
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_BODY-N_ELEM",ninp,nout) ) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) mass
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MASS",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) Ic_imposed
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC_IMPOSED",ninp,nout)) return
       if (Ic_imposed == 1) then
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Ic(1,1),Ic(1,2),Ic(1,3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(1,1-3)",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Ic(2,1),Ic(2,2),Ic(2,3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(2,1-3)",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Ic(3,1),Ic(3,2),Ic(3,3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"IC(3,1-3)",ninp,nout)) return
       endif
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) alfa(1),alfa(2),alfa(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALFA",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) x_rotC(1),x_rotC(2),x_rotC(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_ROTC",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) u_CM(1),u_CM(2),u_CM(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"U_CM",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) omega(1),omega(2),omega(3)
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"OMEGA",ninp,nout)) return
       call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
       read (ainp,*,iostat=ioerr) imposed_kinematics,n_records
       if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS",ninp,nout)) return       
!Assignation to the body parameters 
       body_arr(Id_body)%n_elem = n_elem
       body_arr(Id_body)%mass = mass
       body_arr(Id_body)%x_CM = x_CM
       body_arr(Id_body)%Ic_imposed = Ic_imposed
       if (Ic_imposed == 1) then
          body_arr(Id_body)%Ic = Ic
          else
             body_arr(Id_body)%Ic = 0.
       endif   
       body_arr(Id_body)%alfa = alfa
       body_arr(Id_body)%x_rotC = x_rotC
       body_arr(Id_body)%u_CM = u_CM
       body_arr(Id_body)%omega = omega
       body_arr(Id_body)%imposed_kinematics = imposed_kinematics
       body_arr(Id_body)%n_records = n_records
!Writing on the log file
       if ( ncord > 0 ) then
          if ( nout > 0 ) then
             write (nout,"(1x,a,1p,i12)")    "body:.......................",Id_body
             write (nout,"(1x,a,1p,e12.4)")  "mass:.......................",mass
             write (nout,"(1x,a,1p,3e12.4)") "x_CM:.......................",x_CM
             write (nout,"(1x,a,1p,i12)")    "IC_imposed:.................",Ic_imposed
             if (Ic_imposed == 1) then
             write (nout,"(1x,a,1p,3e12.4)") "Ic(1,1-3):..................",Ic(1,1),Ic(1,2),Ic(1,3)
             write (nout,"(1x,a,1p,3e12.4)") "Ic(2,1-3):..................",Ic(2,1),Ic(2,2),Ic(2,3)
             write (nout,"(1x,a,1p,3e12.4)") "Ic(3,1-3):..................",Ic(3,1),Ic(3,2),Ic(3,3)
             endif
!AA504             
             write (nout,"(1x,a,1p,3e12.4)") "alfa:.......................",alfa(1),alfa(2),alfa(3) 
             write (nout,"(1x,a,1p,3e12.4)") "x_rotC:.....................",x_rotC 
             write (nout,"(1x,a,1p,3e12.4)") "u_CM:.......................",u_CM
             write (nout,"(1x,a,1p,3e12.4)") "omega:......................",omega
             write (nout,"(1x,a,1p,i12)")    "imposed_kinematics:.........",imposed_kinematics    
             write (nout,"(1x,a,1p,i12)")    "n_records:..................",n_records
             write (nout,"(1x,a)")  " "
          end if
       endif
!Allocating body elements
       if (ncord>0) then
          else
          allocate(body_arr(Id_body)%elem(body_arr(Id_body)%n_elem))
          allocate(body_arr(Id_body)%body_kinematics(n_records,7))
       endif
!Reading the eventual imposed kinematics
       do j=1,n_records
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) body_arr(Id_body)%body_kinematics(j,:) 
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY_KINEMATICS_RECORDS",ninp,nout)) return            
       enddo       
! Reading element parameters                     
       do j=1,n_elem
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) Id_elem
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ID_ELEM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) L_geom(1),L_geom(2),L_geom(3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"L_GEOM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) x_CM(1),x_CM(2),x_CM(3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"X_CM_ELEM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) alfa(1),alfa(2),alfa(3)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALFA_ELEM",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) normal_act(1),normal_act(2),normal_act(3),normal_act(4),normal_act(5),normal_act(6)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"NORMAL_ACT",ninp,nout)) return
          call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
          read (ainp,*,iostat=ioerr) mass_deact(1),mass_deact(2),mass_deact(3),mass_deact(4),mass_deact(5),mass_deact(6)
          if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MASS_DEACT",ninp,nout)) return
!Assignation of the element parameters
          body_arr(Id_body)%elem(Id_elem)%L_geom = L_geom
          body_arr(Id_body)%elem(Id_elem)%x_CM = x_CM
          body_arr(Id_body)%elem(Id_elem)%alfa = alfa
          body_arr(Id_body)%elem(Id_elem)%normal_act = normal_act
          body_arr(Id_body)%elem(Id_elem)%mass_deact = mass_deact          
!Writing on the log file
          if ( ncord > 0 ) then
             if ( nout > 0 ) then
                write (nout,"(1x,a,1p,i12)")    "element:....................",Id_elem
                write (nout,"(1x,a,1p,3e12.4)") "L_geom:.....................",L_geom
                write (nout,"(1x,a,1p,3e12.4)") "x_CM_elem:..................",x_CM
!AA504                
                write (nout,"(1x,a,1p,3e12.4)") "alfa_elem:..................",alfa(1),alfa(2),alfa(3)    
                write (nout,"(1x,a,1p,6i12)")   "normal_act:.................",normal_act    
                write (nout,"(1x,a,1p,6e12.4)") "mass_deact:.................",mass_deact                                          
                write (nout,"(1x,a)")  " "
             end if
          endif          
       end do
    enddo         
    call ReadRiga (ainp,comment,nrighe,ioerr,ninp)
    if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"BODY DYNAMICS DATA",ninp,nout)) return
 end do

return
end subroutine ReadBodyDynamics
!---split

