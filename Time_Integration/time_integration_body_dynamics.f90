!cfile time_integration_body_dynamics.f90
!AA501b the whole subroutine 
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : time_integration_body_dynamics
!
! Creation      : 12Jul12 (Amicarelli-Agate)
!
! Last updating : September 20, 2012
!
!************************************************************************************
! Module purpose : Euler time integration for body dynamics
!
! Calling routines: Loop_Irre_2D,Loop_Irre_3D
!
! Called routines: /
!
!************************************************************************************

  subroutine time_integration_body_dynamics(dtvel)  

! Used modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations 
  implicit none
  double precision, intent(in) :: dtvel
  integer(4)       :: i,npi,j
  double precision :: mod_normal,aux_umax
  double precision :: vec_temp(3),vec2_temp(3),vec3_temp(3),domega_dt(3),aux_vel(3)
  double precision,allocatable,dimension(:,:) :: teta

!Allocations
  allocate(teta(n_bodies,3))

! Statements 

! Loop over bodies (body dynamics)
!$omp parallel do default(none) private(i,vec_temp,vec2_temp,vec3_temp,domega_dt,aux_umax) &
!$omp shared(n_bodies,body_arr,dt,ncord,teta,dtvel,tempo)
  do i=1,n_bodies
!staggered parameters (velocity and angular velocity)

     if (body_arr(i)%imposed_kinematics == 0 ) then
! Computed kinematics       
        if (ncord == 2) then
           body_arr(i)%Force(2) = zero
           body_arr(i)%Moment(1) = zero
           body_arr(i)%Moment(3) = zero
        endif
        body_arr(i)%u_CM(:) = body_arr(i)%u_CM(:) + body_arr(i)%Force(:) / body_arr(i)%mass * dtvel 
        if (ncord == 3) then
           vec_temp(1) = dot_product(body_arr(i)%Ic(1,:),body_arr(i)%omega)
           vec_temp(2) = dot_product(body_arr(i)%Ic(2,:),body_arr(i)%omega)
           vec_temp(3) = dot_product(body_arr(i)%Ic(3,:),body_arr(i)%omega)
           call Vector_Product(body_arr(i)%omega,vec_temp,vec3_temp,3)
           vec2_temp = body_arr(i)%Moment(:) - vec3_temp(:)
        endif
        if (ncord == 2) vec2_temp = body_arr(i)%Moment(:) 
        domega_dt(1) = dot_product(body_arr(i)%Ic_inv(1,:),vec2_temp)
        domega_dt(2) = dot_product(body_arr(i)%Ic_inv(2,:),vec2_temp)
        domega_dt(3) = dot_product(body_arr(i)%Ic_inv(3,:),vec2_temp)
        body_arr(i)%omega(:) = body_arr(i)%omega(:) + domega_dt(:) * dtvel
        if (ncord == 2) then
           body_arr(i)%x_CM(2) = zero 
           body_arr(i)%u_CM(2) = zero
           body_arr(i)%omega(1) = zero
           body_arr(i)%omega(3) = zero
        endif
        else
! Imposed kinematics
           do j=1,body_arr(i)%n_records
              if (body_arr(i)%body_kinematics(j,1)>=tempo) then
                 if (body_arr(i)%body_kinematics(j,1)==tempo) then
                    body_arr(i)%u_CM(:) = body_arr(i)%body_kinematics(j,2:4)
                    body_arr(i)%omega(:) = body_arr(i)%body_kinematics(j,5:7)
                    else
                    body_arr(i)%u_CM(:) = body_arr(i)%body_kinematics(j-1,2:4) + &
                                         (body_arr(i)%body_kinematics(j,2:4)-body_arr(i)%body_kinematics(j-1,2:4))/ &
                                         (body_arr(i)%body_kinematics(j,1)-body_arr(i)%body_kinematics(j-1,1)) * &
                                         (tempo-body_arr(i)%body_kinematics(j-1,1))
                    body_arr(i)%omega(:) = body_arr(i)%body_kinematics(j-1,5:7) + &
                                          (body_arr(i)%body_kinematics(j,5:7)-body_arr(i)%body_kinematics(j-1,5:7))/ &
                                         (body_arr(i)%body_kinematics(j,1)-body_arr(i)%body_kinematics(j-1,1)) * &
                                         (tempo-body_arr(i)%body_kinematics(j-1,1))                                      
                 endif
                 exit
              endif
           enddo   
     endif
!     
     
!non-staggered parameters     
     body_arr(i)%x_CM(:) = body_arr(i)%x_CM(:) + body_arr(i)%u_CM(:) * dt 
     teta(i,:) =  body_arr(i)%omega(:) * dt 
     body_arr(i)%alfa(:) = body_arr(i)%alfa(:) + teta(i,:) 
! Initializing umax (maximum particle velocity of the body)
     body_arr(i)%umax = zero 
  end do
!$omp end parallel do

! Loop over body particles (static kinematics)
!$omp parallel do default(none) private(npi,vec_temp,mod_normal,vec2_temp,aux_vel) &
!$omp shared(n_body_part,body_arr,bp_arr,dt,ncord,teta,dtvel)
  do npi=1,n_body_part
! staggered parameter  
     call Vector_Product(body_arr(bp_arr(npi)%body)%omega,bp_arr(npi)%rel_pos,vec_temp,3)
     aux_vel(:) = bp_arr(npi)%vel(:) 
     bp_arr(npi)%vel(:) = body_arr(bp_arr(npi)%body)%u_CM(:) + vec_temp(:)
     if (ncord == 2) bp_arr(npi)%vel(2) = zero 
     bp_arr(npi)%acc(:) = (bp_arr(npi)%vel(:)-aux_vel(:))/dtvel
! non-staggered parameters     
     vec2_temp(:) = teta(bp_arr(npi)%body,:)
     call vector_rotation(bp_arr(npi)%rel_pos,vec2_temp)
     if (ncord == 2) bp_arr(npi)%rel_pos(2) = zero     
     bp_arr(npi)%pos(:) = bp_arr(npi)%rel_pos(:) + body_arr(bp_arr(npi)%body)%x_CM(:)
     call vector_rotation(bp_arr(npi)%normal,vec2_temp)
     mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
     if (mod_normal > one) bp_arr(npi)%normal(:) = bp_arr(npi)%normal(:) / mod_normal    
     if (ncord == 2) then
        bp_arr(npi)%rel_pos(2) = zero
        bp_arr(npi)%normal(2) = zero
     endif 
  end do
!$omp end parallel do
!

!Updating max velocity within every body
  do npi=1,n_body_part
     aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
     body_arr(bp_arr(npi)%body)%umax = max(body_arr(bp_arr(npi)%body)%umax,aux_umax)   
  end do

! Part good for RK1, not for RK1stag
!! Loop over body particles (kinematics)
!!$omp parallel do default(none) private(npi) shared(n_body_part,body_arr,bp_arr,dt,ncord)
!  do npi=1,n_body_part
!     bp_arr(npi)%pos(:) = bp_arr(npi)%pos(:) + bp_arr(npi)%vel(:) * dt
!     if (ncord == 2) bp_arr(npi)%pos(2) = zero 
!  end do
!!$omp end parallel do

!! Loop over bodies (body dynamics)
!!$omp parallel do default(none) private(i,vec_temp,vec2_temp,vec3_temp,domega_dt,aux_umax) shared(n_bodies,body_arr,dt,ncord,teta)
!  do i=1,n_bodies
!     if (ncord == 2) then
!        body_arr(i)%Force(2) = zero
!        body_arr(i)%Moment(1) = zero
!        body_arr(i)%Moment(3) = zero
!     endif
!     body_arr(i)%x_CM(:) = body_arr(i)%x_CM(:) + body_arr(i)%u_CM(:) * dt 
!     teta(i,:) =  body_arr(i)%omega(:) * dt 
!     body_arr(i)%alfa(:) = body_arr(i)%alfa(:) + teta(i,:) 
!     body_arr(i)%u_CM(:) = body_arr(i)%u_CM(:) + body_arr(i)%Force(:) / body_arr(i)%mass * dt 
!     if (ncord == 3) then
!        vec_temp(1) = dot_product(body_arr(i)%Ic(1,:),body_arr(i)%omega)
!        vec_temp(2) = dot_product(body_arr(i)%Ic(2,:),body_arr(i)%omega)
!        vec_temp(3) = dot_product(body_arr(i)%Ic(3,:),body_arr(i)%omega)
!        call Vector_Product(body_arr(i)%omega,vec_temp,vec3_temp,3)
!        vec2_temp = body_arr(i)%Moment(:) - vec3_temp(:)
!     endif
!     if (ncord == 2) vec2_temp = body_arr(i)%Moment(:) 
!     domega_dt(1) = dot_product(body_arr(i)%Ic_inv(1,:),vec2_temp)
!     domega_dt(2) = dot_product(body_arr(i)%Ic_inv(2,:),vec2_temp)
!     domega_dt(3) = dot_product(body_arr(i)%Ic_inv(3,:),vec2_temp)
!     body_arr(i)%omega(:) = body_arr(i)%omega(:) + domega_dt(:) * dt
!     if (ncord == 2) then
!        body_arr(i)%x_CM(2) = zero 
!        body_arr(i)%u_CM(2) = zero
!        body_arr(i)%omega(1) = zero
!        body_arr(i)%omega(3) = zero
!     endif
!! Initializing umax (maximum particle velocity of the body)
!     body_arr(i)%umax = zero 
!  end do
!!$omp end parallel do
!
!! Loop over body particles (static kinematics)
!!$omp parallel do default(none) private(npi,vec_temp,mod_normal,vec2_temp,aux_umax) shared(n_body_part,body_arr,bp_arr,dt,ncord,teta)
!  do npi=1,n_body_part
!     call Vector_Product(body_arr(bp_arr(npi)%body)%omega,bp_arr(npi)%rel_pos,vec_temp,3)
!     bp_arr(npi)%vel(:) = body_arr(bp_arr(npi)%body)%u_CM(:) + vec_temp(:)
!     bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) - body_arr(bp_arr(npi)%body)%x_CM(:)
!     vec2_temp(:) = teta(bp_arr(npi)%body,:)
!     call vector_rotation(bp_arr(npi)%normal,vec2_temp)
!     mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
!     if (mod_normal > 1.) bp_arr(npi)%normal(:) = bp_arr(npi)%normal(:) / mod_normal    
!     if (ncord == 2) then
!        bp_arr(npi)%vel(2) = zero 
!        bp_arr(npi)%rel_pos(2) = zero
!        bp_arr(npi)%normal(2) = zero 
!     endif 
!! Update of umax
!     aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
!     body_arr(bp_arr(npi)%body)%umax = max(body_arr(bp_arr(npi)%body)%umax,aux_umax)   
!  end do
!!$omp end parallel do

!Deallocations
  deallocate(teta)

  return
  end subroutine time_integration_body_dynamics
!---split

