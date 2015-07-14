!cfile Body_dynamics.f90
!AA501b all the subroutines
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: RHS_body_dynamics
!
!AA504 sub start
! Versions: 
! 01  Amicarelli      13nov12       (creation) RHS equations for Body dynamics (body transport scheme)
! 02  Amicarelli      14dic12       Body dynamics loop correction 
! 03  Amicarelli      08apr14       (v5.04) Minor corrections/improvements for: inter-body impingements (particle impact velocities instead of body impact velocities);
!                                   printing close impacts; gravity flag. v504.      
!AA504 end
!
!************************************************************************************
! Module purpose : to estimate the Right Hand Sides of the balance equations for 
!                  body dynamics
!
! Calling routines: Loop_Irre_3D,Loop_Irre_2D
!
! Called subroutines: Vector_Product,FindCloseBoundarySides2D,FindCloseBoundaryFaces3D,
!                     Matrix_Inversion_3x3,Matrix_Inversion_2x2
!
!************************************************************************************
!AA601 changed the name of this subroutine everywhere
  subroutine RHS_body_dynamics

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
!AA504 sub  
  integer(4) :: npartint,i,j,npi,npj,Ncb,Nfzn,aux,nbi,npk,k,nbj,nbk,n_interactions,aux2,aux3,test,aux_locx_min,aux_locx_max,aux_int
  double precision :: c2,k_masses,r_per,r_par,temp_dden,temp_acc,alfa_boun,aux_impact_vel,aux4,pres_mir
  double precision :: f_pres(3),temp(3),r_par_vec(3),f_coll_bp_bp(3),f_coll_bp_boun(3),dvar(3),pos_aux(3),normal_plane(3)
  double precision :: u_rel(3),x_rel(3),aux_acc(3),aux_vec(3),aux_vec2(3),loc_pos(3),aux_locx_vert(3)
  double precision :: aux_mat(3,3)
  double precision,dimension(:,:,:),allocatable :: Force,Moment
  double precision,dimension(:,:),allocatable :: Force_mag_sum,r_per_min,aux_gravity
!AA504
  character(255) :: file_name_test
    
! External functions
  double precision, external :: Gamma_boun

!Allocations
  if (ncord==2) aux = n_bodies+NumBSides
  if (ncord==3) aux = n_bodies+NumFacce
  allocate(Force(n_bodies,aux,3))
  allocate(Moment(n_bodies,aux,3))
  allocate(Force_mag_sum(n_bodies,aux))
  allocate(r_per_min(n_bodies,aux))
  allocate(aux_gravity(n_bodies,3))
  
!Initializations
  Force = 0.
  Moment = 0.
  Force_mag_sum = 0.
  r_per_min = 1000000.
  aux2 = 0

! Statements

! Updating pressure of the body particles
  call body_pressure_mirror

! Contributions to fluid dynamics momentum (discretized semi-analytic approach: a mirror particle technique)
! Loop over body particles (Not in parallel because of critical sections)
  do npi=1,n_body_part
     do j=1,nPartIntorno_bp_f(npi)
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno_bp_f(npartint)
        temp_acc =  (pg(npj)%pres+bp_arr(npi)%pres)/(pg(npj)%dens*pg(npj)%dens)
        pg(npj)%acc(:) = pg(npj)%acc(:) + ( - pg(npj)%mass/(dx_dxbodies**ncord) * temp_acc * &
                         (-rag_bp_f(:,npartint)) * KerDer_bp_f_Gal(npartint) )    
     end do
  end do

! Loop over the transported bodies (gravity forces and Ic initialization)
!$omp parallel do default(none) private(i) shared(n_bodies,body_arr,Domain,it_start,it_corrente,ncord,aux_gravity)
  do i=1,n_bodies
  
! No computations for imposed kinematics
     if (body_arr(i)%imposed_kinematics == 0 ) then
     
     aux_gravity(i,:) = body_arr(i)%mass * Domain%grav(:)
     body_arr(i)%Force = 0.
     body_arr(i)%Moment = 0.
     if ((body_arr(i)%Ic_imposed == 0) .and. ((ncord == 3).or.((it_start+1) == it_corrente)) ) then
        body_arr(i)%Ic = 0.
     endif
     
     endif ! No computations for imposed kinematics
     
  end do
!$omp end parallel do

! Loop over the body particles  
!omp parallel do default(none) private(npi,f_pres,temp,j,npartint,npj,temp_acc,dvar,temp_dden,c2,r_per,r_par_vec) &
!omp private(r_par,k_masses,f_coll_bp_bp,f_coll_bp_boun,pos_aux,normal_plane,aux2,u_rel,x_rel,aux_impact_vel,k,nbi,i,aux_vec) &
!omp private(aux_vec2,aux4,aux_mat,loc_pos,test,aux_locx_vert,aux_locx_min,aux_locx_max) &
!omp shared(n_body_part,bp_arr,body_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f,pg,KerDer_bp_f_Gal,rag_bp_f,dx_dxbodies) &
!omp shared(nPartIntorno_bp_bp,PartIntorno_bp_bp,Med,rag_bp_bp,Domain,ncord,NumFacce,BoundaryFace,BoundarySide,NumBSides) &
!omp shared(impact_vel,Vertice,n_surf_body_part,surf_body_part,Force,imping_body_grav,aux_gravity,Moment,Force_mag_sum) &
!omp shared(r_per_min,Tratto,n_bodies)
  do npi=1,n_body_part

! No computations for imposed kinematics
     if (body_arr(bp_arr(npi)%body)%imposed_kinematics == 0 ) then

! Auxiliary variable to indicate surface body particles     
     if (bp_arr(npi)%area > 0.) then ! only body particles at the body surface 

!Computation of the ID of the surface body particles
!        i = 0
!        aux2 = 0
!        do while (i<npi) 
!           i = i+1 
!           if (bp_arr(i)%area > 0.) aux2 = aux2+1
!        enddo
        aux2 = aux2+1 
     
! Fluid pressure forces
        f_pres(:) = bp_arr(npi)%pres * bp_arr(npi)%area * bp_arr(npi)%normal(:)
        body_arr(bp_arr(npi)%body)%Force(:) = body_arr(bp_arr(npi)%body)%Force(:) + f_pres(:)
        call Vector_Product(bp_arr(npi)%rel_pos(:),f_pres(:),temp(:),3)
        body_arr(bp_arr(npi)%body)%Moment(:) = body_arr(bp_arr(npi)%body)%Moment(:) + temp(:)

! Loop over the neighbouring body particles, belonging to other bodies (inter-body impacts, partial contributions) 
        do j=1,nPartIntorno_bp_bp(aux2)
           npartint = (aux2-1)* NMAXPARTJ + j
           npj = PartIntorno_bp_bp(npartint)
           r_per = abs(rag_bp_bp(1,npartint)*bp_arr(npj)%normal(1)+ &
                       rag_bp_bp(2,npartint)*bp_arr(npj)%normal(2)+ &
                       rag_bp_bp(3,npartint)*bp_arr(npj)%normal(3)) 
           r_par_vec(:) = - rag_bp_bp(:,npartint) - r_per * bp_arr(npj)%normal(:)  ! (note: rag=-r)
           r_par = dsqrt(r_par_vec(1)*r_par_vec(1)+r_par_vec(2)*r_par_vec(2)+r_par_vec(3)*r_par_vec(3))
           if ( (r_per>0.) .and. ((r_par/(Domain%dd/dx_dxbodies)) <= 1.) ) then
!AA504 sub start              
               u_rel(:) = bp_arr(npi)%vel(:) - bp_arr(npj)%vel(:) 
               x_rel(:) = - rag_bp_bp(:,npartint) / dsqrt(dot_product(rag_bp_bp(:,npartint),rag_bp_bp(:,npartint)))
!AA504 sub end              
              aux_impact_vel = dot_product(u_rel,x_rel)
              if ( (impact_vel(aux2,bp_arr(npj)%body)==0.) .or. (impact_vel(aux2,bp_arr(npj)%body)<aux_impact_vel) ) then
                 k = 0
                 nbi = 1
                 do while (nbi.ne.bp_arr(npi)%body) 
                       k = k + body_arr(nbi)%npart
                       nbi = nbi +1
                 enddo
                 do npk=(k+1),(k+body_arr(bp_arr(npi)%body)%npart)
!AA502
!                   do aux3=i,n_surf_body_part
                    do aux3=1,n_surf_body_part
                       if (surf_body_part(aux3)==npk) then
                          impact_vel(aux3,bp_arr(npj)%body) = body_arr(bp_arr(npi)%body)%umax+body_arr(bp_arr(npj)%body)%umax
                          cycle
                       endif
                    enddo
                 enddo
              endif  
              if (impact_vel(aux2,bp_arr(npj)%body)>0.) then 
                 k_masses = (bp_arr(npi)%mass * bp_arr(npj)%mass) / (bp_arr(npi)%mass + bp_arr(npj)%mass)
                 f_coll_bp_bp(:) = - (2*(impact_vel(aux2,bp_arr(npj)%body)**2)/r_per) * k_masses * Gamma_boun(r_per,Domain%h) &
                                   * (1-r_par/(Domain%dd/dx_dxbodies)) * bp_arr(npj)%normal(:)

!AA504 test start
                 if (r_per<((Domain%dd/dx_dxbodies)/10.)) then   
                    write(file_name_test,"(a,a,i8.8,a,i8.8,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_close_impact_',it_corrente,'_',npi,'_',npj,".txt"
                    open (60,file=file_name_test,status="unknown",form="formatted")
                    write (60,'((10x,a),4(11x,a),(2x,a),(9x,a),(6x,a),(4x,a),(9x,a),3(11x,a))') &
                       " time"," npi"," npj"," nbi"," nbj"," imp_vel(m/s)"," r_per"," k_masses"," Gamma_boun"," r_par", &
                       " n_x"," n_y"," n_z"
                    write (60,'((g14.7,1x),4(i14,1x),8(g14.7,1x))') tempo,npi,npj,bp_arr(npi)%body,bp_arr(npj)%body,impact_vel(aux2,bp_arr(npj)%body),r_per,k_masses, &
                       Gamma_boun(r_per,Domain%h),r_par,bp_arr(npj)%normal(1),bp_arr(npj)%normal(2),bp_arr(npj)%normal(3)
                   close(60)
                 endif
!AA504 test end
                 
                 Force(bp_arr(npi)%body,bp_arr(npj)%body,:) = Force(bp_arr(npi)%body,bp_arr(npj)%body,:) + f_coll_bp_bp(:)
! zeroing gravity component perpendicular to the normal    
                 if (imping_body_grav==0) then
!AA504 sub start  
                     aux_gravity(bp_arr(npi)%body,:) = 0.d0
! interesting test
!                    aux_vec(:) = dot_product(aux_gravity(bp_arr(npi)%body,:),bp_arr(npj)%normal) * bp_arr(npj)%normal(:)
!                    aux_gravity(bp_arr(npi)%body,:) = aux_gravity(bp_arr(npi)%body,:) - aux_vec(:)
!AA504 sub end
                 endif
!              
                 call Vector_Product(bp_arr(npi)%rel_pos,f_coll_bp_bp,temp,3)
                 Moment(bp_arr(npi)%body,bp_arr(npj)%body,:) = Moment(bp_arr(npi)%body,bp_arr(npj)%body,:) + temp(:)
                 Force_mag_sum(bp_arr(npi)%body,bp_arr(npj)%body) = Force_mag_sum(bp_arr(npi)%body,bp_arr(npj)%body) + &
                                     (1/r_per) * k_masses * Gamma_boun(r_per,Domain%h) * (1-r_par/(Domain%dd/dx_dxbodies))
                 r_per_min(bp_arr(npi)%body,bp_arr(npj)%body) = min(r_per,r_per_min(bp_arr(npi)%body,bp_arr(npj)%body)) 
              endif
           endif                            
        end do

! Loop over boundaries (body-boundary impacts, partial contributions)
!3D case 
        if (ncord == 3) then
           do j=1,NumFacce
           
!AA501btest           
!              if (BoundaryFace(j)%stretch == 1) then
              if ((Tratto(BoundaryFace(j)%stretch)%tipo == "fixe") .or. (Tratto(BoundaryFace(j)%stretch)%tipo == "tapi")) then

!AA501btest
                 aux_vec2(:) = bp_arr(npi)%pos(:) - BoundaryFace(j)%Node(1)%GX(:)
                 aux4 = dot_product(BoundaryFace(j)%T(:,3),aux_vec2)
                 if (aux4>=0.) then
              
                    call dis_point_plane(bp_arr(npi)%pos,BoundaryFace(j)%Node(1)%GX,BoundaryFace(j)%Node(2)%GX, &
                                         BoundaryFace(j)%Node(3)%GX,r_per,normal_plane)
                                      
!AA501btest start
!                 if ( (r_per>0.) .and. (r_per <= (2.*Domain%h)) ) then
                    aux_mat(:,1) = BoundaryFace(j)%T(:,1)
                    aux_mat(:,2) = BoundaryFace(j)%T(:,2)
                    aux_mat(:,3) = BoundaryFace(j)%T(:,3)
                    call reference_system_change(bp_arr(npi)%pos,BoundaryFace(j)%Node(4)%GX,aux_mat,loc_pos)
                    call point_inout_polygone(loc_pos(1:2),BoundaryFace(j)%nodes,BoundaryFace(j)%Node(1)%LX(1:2), &
                            BoundaryFace(j)%Node(2)%LX(1:2),BoundaryFace(j)%Node(3)%LX(1:2), &
                            BoundaryFace(j)%Node(4)%LX(1:2),test)
                    if ( (r_per>0.) .and. (r_per <= (2.*Domain%h)) .and. (test==1) ) then     
!AA501btest end                                     

                       aux_impact_vel = -(dot_product(normal_plane,bp_arr(npi)%vel))
                       if ( (impact_vel(aux2,n_bodies+j)==0.) .or. &
                            (impact_vel(aux2,n_bodies+j)<aux_impact_vel) ) then
                             impact_vel(aux2,n_bodies+j) = aux_impact_vel
                       endif
                       if (impact_vel(aux2,n_bodies+j)>0.) then 
                          k_masses = bp_arr(npi)%mass 
                          f_coll_bp_boun(:) = (2.*(impact_vel(aux2,n_bodies+j)**2)/r_per) * k_masses * Gamma_boun(r_per,Domain%h) &
                                              * normal_plane(:)
                          Force(bp_arr(npi)%body,n_bodies+j,:) = Force(bp_arr(npi)%body,n_bodies+j,:) + f_coll_bp_boun(:) 
! zeroing gravity component perpendicular to the normal
                          if (imping_body_grav==0) then  
!AA504 sub start
                              aux_gravity(bp_arr(npi)%body,:) = 0.d0
! interesting test            
!                             aux_vec(:) = dot_product(aux_gravity(bp_arr(npi)%body,:),normal_plane) * normal_plane(:)
!                             aux_gravity(bp_arr(npi)%body,:) = aux_gravity(bp_arr(npi)%body,:) - aux_vec(:)
!AA504 sub end
                          endif
!                       
                          call Vector_Product(bp_arr(npi)%rel_pos,f_coll_bp_boun,temp,3)
                          Moment(bp_arr(npi)%body,n_bodies+j,:) = Moment(bp_arr(npi)%body,n_bodies+j,:) + temp(:)
                          Force_mag_sum(bp_arr(npi)%body,n_bodies+j) = Force_mag_sum(bp_arr(npi)%body,n_bodies+j) + &
                                                                       (1/r_per) * k_masses * Gamma_boun(r_per,Domain%h) 
                          r_per_min(bp_arr(npi)%body,n_bodies+j) = min(r_per,r_per_min(bp_arr(npi)%body,n_bodies+j))
                          if (Gamma_boun(r_per,Domain%h)<=0.) impact_vel(aux2,n_bodies+j) = 0.
                       endif
                    endif

!AA501btest
                 endif  

              endif
           end do 
        endif
!2D case
        if (ncord == 2) then
           do j=1,NumBSides
           
!AA501btest           
!              if (BoundarySide(j)%stretch == 1) then

              if ((BoundarySide(j)%tipo == "fixe") .or. (BoundarySide(j)%tipo == "tapi")) then  
                 aux_vec2(:) = bp_arr(npi)%pos(:) - Vertice(1,BoundarySide(j)%Vertex(1))
                 aux4 = dot_product(BoundarySide(j)%T(:,3),aux_vec2)
                 if (aux4>=0.) then
                    pos_aux(1) =  Vertice(1,BoundarySide(j)%Vertex(1))
                    pos_aux(2) = 0.
                    pos_aux(3) =  Vertice(3,BoundarySide(j)%Vertex(1))
                    call dis_point_plane(bp_arr(npi)%pos,Vertice(:,BoundarySide(j)%Vertex(1)), &
                                         Vertice(:,BoundarySide(j)%Vertex(2)),pos_aux,r_per,normal_plane)
                                      
!AA501btest start
!                 if ( (r_per>0.) .and. (r_per <= (2.*Domain%h)) ) then
                    aux_mat(:,1) = BoundarySide(j)%T(:,1)
                    aux_mat(2,1) = 0.
                    aux_mat(2,2) = 1.
                    aux_mat(2,3) = 0.
                    aux_mat(:,3) = BoundarySide(j)%T(:,3) 
                    call reference_system_change(bp_arr(npi)%pos,Vertice(:,BoundarySide(j)%Vertex(1)),aux_mat,loc_pos)
                    call reference_system_change(Vertice(:,BoundarySide(j)%Vertex(2)),Vertice(:,BoundarySide(j)%Vertex(1)), &
                                                 aux_mat,aux_locx_vert)
                    aux_locx_min = min(zero,aux_locx_vert(1))
                    aux_locx_max = max(zero,aux_locx_vert(1))
                    if ((loc_pos(1)>=aux_locx_min) .and. &
                        (loc_pos(1)<=aux_locx_max)) then
                       test = 1
                       else
                          test = 0
                    endif
                    if ( (r_per>0.) .and. (r_per <= (2.*Domain%h)) .and. (test==1) ) then     
!AA501btest end                                       
                                      
                       aux_impact_vel = -(dot_product(normal_plane,bp_arr(npi)%vel))
                       if ( (impact_vel(aux2,n_bodies+j)==0.) .or. &
                            (impact_vel(aux2,n_bodies+j)<aux_impact_vel) ) then
                             impact_vel(aux2,n_bodies+j) = aux_impact_vel
                       endif
                       if (impact_vel(aux2,n_bodies+j)>0.) then 
                          k_masses = bp_arr(npi)%mass 
                          f_coll_bp_boun(:) = (2.*(impact_vel(aux2,n_bodies+j)**2)/r_per) * k_masses * Gamma_boun(r_per,Domain%h) &
                                              * normal_plane(:)
                          Force(bp_arr(npi)%body,n_bodies+j,:) = Force(bp_arr(npi)%body,n_bodies+j,:) + f_coll_bp_boun(:) 
! zeroing gravity component perpendicular to the normal    
                          if (imping_body_grav==0) then
!AA504 sub start  
                              aux_gravity(bp_arr(npi)%body,:) = 0.d0
! interesting test
!                             aux_vec(:) = dot_product(aux_gravity(bp_arr(npi)%body,:),normal_plane) * normal_plane(:)
!                             aux_gravity(bp_arr(npi)%body,:) = aux_gravity(bp_arr(npi)%body,:) - aux_vec(:)  
!AA504 sub end
                          endif                    
!
                          call Vector_Product(bp_arr(npi)%rel_pos,f_coll_bp_boun,temp,3)
                          Moment(bp_arr(npi)%body,n_bodies+j,:) = Moment(bp_arr(npi)%body,n_bodies+j,:) + temp(:)
                          Force_mag_sum(bp_arr(npi)%body,n_bodies+j) = Force_mag_sum(bp_arr(npi)%body,n_bodies+j) + &
                                                                       (1/r_per) * k_masses * Gamma_boun(r_per,Domain%h) 
                          r_per_min(bp_arr(npi)%body,n_bodies+j) = min(r_per,r_per_min(bp_arr(npi)%body,n_bodies+j))
                          if (Gamma_boun(r_per,Domain%h)<=0.) impact_vel(aux2,n_bodies+j) = 0.
                       endif
                    endif
                 endif
              endif
           end do 
        endif
!        
     endif ! only body particle at the body surface
     
     endif ! No computations for imposed kinematics
     
  end do

! Contributions to the moment of inertia
 do npi=1,n_body_part
 
! No computations for imposed kinematics
    if (body_arr(bp_arr(npi)%body)%imposed_kinematics == 0 ) then
    
    if (body_arr(bp_arr(npi)%body)%Ic_imposed == 0) then
! Contributions to the moment of inertia coefficients (diagonal elements)
    if (ncord == 3) then
       body_arr(bp_arr(npi)%body)%Ic(1,1) = body_arr(bp_arr(npi)%body)%Ic(1,1) + bp_arr(npi)%mass * &
                                            (bp_arr(npi)%rel_pos(2)**2 + bp_arr(npi)%rel_pos(3)**2)
       body_arr(bp_arr(npi)%body)%Ic(2,2) = body_arr(bp_arr(npi)%body)%Ic(2,2) + bp_arr(npi)%mass * &
                                            (bp_arr(npi)%rel_pos(1)**2 + bp_arr(npi)%rel_pos(3)**2)
       body_arr(bp_arr(npi)%body)%Ic(3,3) = body_arr(bp_arr(npi)%body)%Ic(3,3) + bp_arr(npi)%mass * &
                                            (bp_arr(npi)%rel_pos(1)**2 + bp_arr(npi)%rel_pos(2)**2)
! Contributions to the products of inertia (off-diagonal elements)
       body_arr(bp_arr(npi)%body)%Ic(1,2) = body_arr(bp_arr(npi)%body)%Ic(1,2) - bp_arr(npi)%mass * &
                                            (bp_arr(npi)%rel_pos(1) * bp_arr(npi)%rel_pos(2))
       body_arr(bp_arr(npi)%body)%Ic(1,3) = body_arr(bp_arr(npi)%body)%Ic(1,3) - bp_arr(npi)%mass * &
                                            (bp_arr(npi)%rel_pos(1) * bp_arr(npi)%rel_pos(3))
       body_arr(bp_arr(npi)%body)%Ic(2,3) = body_arr(bp_arr(npi)%body)%Ic(2,3) - bp_arr(npi)%mass * &
                                            (bp_arr(npi)%rel_pos(2) * bp_arr(npi)%rel_pos(3)) 
    endif 
    if ((ncord == 2).and.((it_start+1) == it_corrente)) then
       body_arr(bp_arr(npi)%body)%Ic(2,2) = body_arr(bp_arr(npi)%body)%Ic(2,2) + bp_arr(npi)%mass * &
                                            (bp_arr(npi)%rel_pos(1)**2 + bp_arr(npi)%rel_pos(3)**2)
    endif
    endif 

    endif ! No computations for imposed kinematics
    
 enddo     

! Loop over the transported bodies (global contributions from body-body and boundary-body impacts;
!                                   computation of Ic and its inverse)
!AA504 sub omp directives
!$omp parallel do default(none) private(i,j,k_masses,alfa_boun,aux_int) &
!$omp shared(n_bodies,body_arr,ncord,it_start,it_corrente,aux,r_per_min,Domain,Force_mag_sum,Force,Moment)
 do i=1,n_bodies

! No computations for imposed kinematics
    if (body_arr(i)%imposed_kinematics == 0 ) then

!Forces and moments
    do j=1,aux
       if (r_per_min(i,j)<1000000.) then
          if (j<=n_bodies) then
             k_masses = (body_arr(i)%mass * body_arr(j)%mass) / (body_arr(i)%mass + body_arr(j)%mass)
             else
             k_masses = body_arr(i)%mass 
          endif
          if (Force_mag_sum(i,j)>0.) then
             alfa_boun = (Gamma_boun(r_per_min(i,j),Domain%h)/r_per_min(i,j)*k_masses)/Force_mag_sum(i,j)
             Force(i,j,:) = Force(i,j,:) * alfa_boun
             body_arr(i)%Force(:) = body_arr(i)%Force(:) + Force(i,j,:)
             Moment(i,j,:) = Moment(i,j,:) * alfa_boun
             body_arr(i)%Moment(:) = body_arr(i)%Moment(:) + Moment(i,j,:)
          endif
       endif
    end do

! Ic and its inverse     
    if (ncord == 3) then
       body_arr(i)%Ic(2,1) = body_arr(i)%Ic(1,2)
       body_arr(i)%Ic(3,1) = body_arr(i)%Ic(1,3)
       body_arr(i)%Ic(3,2) = body_arr(i)%Ic(2,3)
       call Matrix_Inversion_3x3(body_arr(i)%Ic,body_arr(i)%Ic_inv,aux_int)
    endif
    if ((ncord == 2).and.((it_start+1) == it_corrente)) then
       body_arr(i)%Ic_inv = 0. 
       body_arr(i)%Ic_inv(2,2) = 1./body_arr(i)%Ic(2,2)
    endif

    endif ! No computations for imposed kinematics
    
  end do
!$omp end parallel do     
  
!Check body-body interactions and eventually zeroing impact velocities; adding gravity 
!$omp parallel do default(none) &
!$omp private(j,k,nbi,npk,nbk,n_interactions,aux2) &
!$omp shared(n_bodies,Moment,Force,body_arr,impact_vel,aux,Domain,n_surf_body_part,surf_body_part,aux_gravity)
  do nbi=1,n_bodies
  
! No computations for imposed kinematics
     if (body_arr(nbi)%imposed_kinematics == 0 ) then  
  
!gravity contribution
     body_arr(nbi)%Force(:) = body_arr(nbi)%Force(:) + aux_gravity(nbi,:) 
!impact velocity  
     n_interactions = 0
     do j=1,aux
        if (nbi==j) cycle
        if ((Moment(nbi,j,1)==0.).and.(Force(nbi,j,1)==0.).and. &
            (Moment(nbi,j,2)==0.).and.(Force(nbi,j,2)==0.).and. &
            (Moment(nbi,j,3)==0.).and.(Force(nbi,j,3)==0.)) then
           if (j<=n_bodies) then
              k = 0
              nbk = 1
              do while (nbk.ne.nbi) 
                    k = k + body_arr(nbk)%npart
                    nbk = nbk +1
              enddo
              do npk=(k+1),(k+body_arr(nbi)%npart)
                 do aux2=1,n_surf_body_part
                    if (surf_body_part(aux2)==npk) then
                       impact_vel(aux2,j) = zero
                       cycle
                    endif
                 enddo
              enddo
           endif
           else
              n_interactions = n_interactions + 1 
        endif
     end do
     
! No computations for imposed kinematics
     endif
     
  enddo
!$omp end parallel do     

!Deallocations
  deallocate(Force)
  deallocate(Moment)
  deallocate(Force_mag_sum)
  deallocate(r_per_min)
  deallocate(aux_gravity)

  return
  end subroutine RHS_body_dynamics
!---split

!cfile body_pressure_mirror.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_pressure_mirror
!
! Creation: Amicarelli-Agate 13nov12
!
!************************************************************************************
! Module purpose : Computation of pressure for body particles (a mirror particle technique)
!
! Calling routines: Loop_Irre_2D,Loop_Irre_3D,RHS_body_dynamics
!
! Called subroutines: /
!
!************************************************************************************
!
  subroutine body_pressure_mirror

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: Sum_W_vol,W_vol,dis,pres_mir,aux
  double precision :: aux_acc(3)

! External functions
  double precision, external :: w

! Statements
! Loop over body particles (discretized semi-analytic approach: a mirror particle technique)
!$omp parallel do default(none) &
!$omp private(npi,Sum_W_vol,j,npartint,npj,aux_acc,pres_mir,dis,W_vol,aux) &
!$omp shared(n_body_part,bp_arr,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f,Domain,pg,rag_bp_f)
  do npi=1,n_body_part
     bp_arr(npi)%pres = 0.
     Sum_W_vol = 0.  
     do j=1,nPartIntorno_bp_f(npi)
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno_bp_f(npartint)

!AA501btest start        
        aux = dsqrt(dot_product(bp_arr(npi)%acc,bp_arr(npi)%acc))
! wall acceleration should be less than 100m/2^2, otherwise an impulse is assumed to occur &
! and the formulation with acc_body is not valid

!AA501btest
        if (aux <= 100.) then
!if (aux <= 10.) then

           aux_acc(:) = Domain%grav(:) - bp_arr(npi)%acc(:)
           else

!AA501btest           
              aux_acc(:) = Domain%grav(:) - (100./aux) * bp_arr(npi)%acc(:)
!   aux_acc(:) = Domain%grav(:) - (10./aux) * bp_arr(npi)%acc(:)

        endif
!AA501btest end        

        pres_mir = pg(npj)%pres - pg(npj)%dens * dot_product(aux_acc,-rag_bp_f(:,npartint))
        dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
        W_vol = w(dis,Domain%h,Domain%coefke) * (pg(npj)%mass/pg(npj)%dens)
        bp_arr(npi)%pres = bp_arr(npi)%pres + pres_mir * W_vol
        Sum_W_vol = Sum_W_vol + W_vol  
                  
     end do
     if (Sum_W_vol > 0.) bp_arr(npi)%pres = bp_arr(npi)%pres / Sum_W_vol      
  end do
 !$omp end parallel do
 
  return
  end subroutine body_pressure_mirror
!---split


!cfile body_particles_to_continuity.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_particles_to_continuity
!
! Creation: Amicarelli-Agate 13nov12
!
!************************************************************************************
! Module purpose : contributions of body particles to the continuity equation
!
! Calling routines: inter_EqCont_3D,inter_EqCont_2D
!
! Called subroutines: /
!
!************************************************************************************

  subroutine body_particles_to_continuity

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!AA501btest
  use files_entities  

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj,k
  double precision :: temp_dden,aux,dis,dis_min,x_min,x_max,y_min,y_max,z_min,z_max,mod_normal,W_vol,sum_W_vol
  double precision :: dvar(3),aux_vec(3),aux_nor(3),aux_vec2(3)
  double precision, external :: w  

! Statements

! Loop over the body particles (maybe in parallel with a crticial section)
!omp parallel do default(none) &
!omp private(npi,j,npartint,npj,dvar,temp_dden) &
!omp shared(n_body_part,nPartIntorno_bp_f,NMAXPARTJ,PartIntorno_bp_f,bp_arr,pg,KerDer_bp_f_cub_spl,rag_bp_f) &
  do npi=1,n_body_part
  
    bp_arr(npi)%vel_mir=0.
    sum_W_vol = 0.
  
! Loop over fluid particles (contributions to fluid particles, discretized semi-analytic approach: mirror particle technique) 
     do j=1,nPartIntorno_bp_f(npi)
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno_bp_f(npartint)
! Continuity equation
! Normal for uSA
        dis_min = 999999999.
!AA501btest start        
        x_min = min(bp_arr(npi)%pos(1),pg(npj)%coord(1))
        x_max = max(bp_arr(npi)%pos(1),pg(npj)%coord(1))
        y_min = min(bp_arr(npi)%pos(2),pg(npj)%coord(2))
        y_max = max(bp_arr(npi)%pos(2),pg(npj)%coord(2))
        z_min = min(bp_arr(npi)%pos(3),pg(npj)%coord(3))
        z_max = max(bp_arr(npi)%pos(3),pg(npj)%coord(3))
!AA501btest end
        aux_nor = 0.
        do k=1,n_body_part
           if (bp_arr(k)%body==bp_arr(npi)%body) then
              aux=dot_product(rag_bp_f(:,npartint),bp_arr(k)%normal)
              if (aux>0.) then
                 if (npi==k) then
                    aux_nor(:) = bp_arr(k)%normal(:)
                    exit
                    else
                       if ( (bp_arr(k)%pos(1)>=x_min) .and. (bp_arr(k)%pos(1)<=x_max) .and. &
                            (bp_arr(k)%pos(2)>=y_min) .and. (bp_arr(k)%pos(2)<=y_max) .and. &
                            (bp_arr(k)%pos(3)>=z_min) .and. (bp_arr(k)%pos(3)<=z_max) ) then
                          call distance_point_line_3D(bp_arr(k)%pos,bp_arr(npi)%pos,pg(npj)%coord,dis)
                          dis_min =min(dis_min,dis)
                          if (dis==dis_min) aux_nor(:) = bp_arr(k)%normal(:)
                       endif
                 endif
              endif
           endif
        end do
        
! In case the normal is not yet defined the normal vector is that of the body particle (at the body surface), 
!    which is closest to the fluid particle
!AA501btest start
        mod_normal = dsqrt(dot_product(aux_nor,aux_nor))
        if (mod_normal==0.) then
           dis_min = 999999999.
           do k=1,n_body_part
              mod_normal=dot_product(bp_arr(k)%normal,bp_arr(k)%normal)
              if (mod_normal>0.) then
                 aux_vec(:) = bp_arr(k)%pos(:) - pg(npj)%coord(:) 
                 dis = dsqrt(dot_product(aux_vec,aux_vec))
                 dis_min =min(dis_min,dis)
                 if (dis==dis_min) aux_nor(:) = bp_arr(k)%normal(:)
              endif
           enddo
        endif

! relative velocity for continuity equation     
        aux_vec2(:) = bp_arr(npi)%vel(:) - pg(npj)%vel(:)
        dvar(:) = aux_nor(:)*2.*dot_product(aux_vec2,aux_nor)
        dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
        W_vol = w(dis,Domain%h,Domain%coefke) * (pg(npj)%mass/pg(npj)%dens)
        bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) + (dvar(:)+pg(npj)%vel(:)) * W_vol
        sum_W_vol = sum_W_vol + W_vol            
       
! Contributions to continuity equations       
        temp_dden = pg(npj)%mass/(dx_dxbodies**ncord) * KerDer_bp_f_cub_spl(npartint) * &
                   ( dvar(1)*(-rag_bp_f(1,npartint)) + dvar(2)*(-rag_bp_f(2,npartint)) + dvar(3)*(-rag_bp_f(3,npartint)) )
        pg(npj)%dden = pg(npj)%dden - temp_dden
     
     end do

     if (sum_W_vol > 0.)   bp_arr(npi)%vel_mir(:) = bp_arr(npi)%vel_mir(:) / sum_W_vol     

  end do
!omp end parallel do
 
  return
  end subroutine body_particles_to_continuity
!---split



!cfile Input_Body_Dynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: Input_Body_Dynamics
!
! Creation: Amicarelli-Agate 16July12
!
!************************************************************************************
! Module purpose : Input management for body dynamics
!
! Calling routines: Gest_Input
!
! Called subroutines: MatrixProduct
!
!AA504 sub start
! Versions: 
! 01  Amicarelli      13nov12       (creation) Input management for the body transport scheme 
! 02  Amicarelli      08apr14       (v5.04) removed error in parallelization (input management)      
!AA504 sub end
!
!************************************************************************************
!
  subroutine Input_Body_Dynamics

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use files_entities

! Declarations
  implicit none
  integer(4) :: i,nbi,npi,j,k,nei,erased_part,ier
  integer(4) :: vec_temp(3)
  double precision :: xmin,xmax,ymin,ymax,zmin,zmax,mod_normal,aux_umax
  double precision :: vec2_temp(3)
  double precision, dimension(:), allocatable  ::  aux_mass
  type (body_particle), dimension(:), allocatable :: aux_bp_arr 

! External functions
  integer(4), external :: ParticleCellNumber

! Allocations
  allocate(aux_mass(n_bodies))

! Initializations
  n_body_part = 0
  aux_mass = 0.
  erased_part = 0
  n_surf_body_part = 0

! Statements
! Loop over the transported bodies
!$omp parallel do default(none) private(i) &
!$omp shared(n_bodies,body_arr,ncord,Domain,n_body_part,dx_dxbodies)
  do i=1,n_bodies
      body_arr(i)%npart = 0
      do j=1,body_arr(i)%n_elem
! Particle length scales for the element
         body_arr(i)%elem(j)%dx(:) = body_arr(i)%elem(j)%L_geom(:) / (int(body_arr(i)%elem(j)%L_geom(:) / (Domain%dd/dx_dxbodies)))
! Number of body particles of the element
         if (ncord == 3) then
            body_arr(i)%elem(j)%npart = int(body_arr(i)%elem(j)%L_geom(1)/body_arr(i)%elem(j)%dx(1)) * &
                                        int(body_arr(i)%elem(j)%L_geom(2)/body_arr(i)%elem(j)%dx(2)) * &
                                        int(body_arr(i)%elem(j)%L_geom(3)/body_arr(i)%elem(j)%dx(3))
         endif
         if (ncord == 2) then
            body_arr(i)%elem(j)%npart = int(body_arr(i)%elem(j)%L_geom(1)/body_arr(i)%elem(j)%dx(1)) * &
                                        int(body_arr(i)%elem(j)%L_geom(3)/body_arr(i)%elem(j)%dx(3))
         endif
! Update of the number of body particles of the body
         body_arr(i)%npart = body_arr(i)%npart + body_arr(i)%elem(j)%npart         
      end do
!AA504 removed line
! Initializing umax (maximum particle velocity of the body)
      body_arr(i)%umax = 0. 
      body_arr(i)%pmax = 0.
   enddo
!$omp end parallel do

!AA504 start
! Incrementing the total number of particles 
  do i=1,n_bodies
     n_body_part = n_body_part + body_arr(i)%npart      
  end do
!AA504 end
  
! Managing body particles  
  allocate (bp_arr(n_body_part))
  npi = 1  

! Loop over the transported bodies (not to be parallelized)
  do nbi=1,n_bodies
  
! Loop over the body elements  
     do nei=1,body_arr(nbi)%n_elem
     vec_temp(:) = int(body_arr(nbi)%elem(nei)%L_geom(:) / body_arr(nbi)%elem(nei)%dx(:))     
     if (ncord == 2) vec_temp(2) = 1     
     do i=1,vec_temp(1)
        do j=1,vec_temp(2)
           do k=1,vec_temp(3)
! Corresponding Body 
              bp_arr(npi)%body = nbi
! Initial relative positions (body reference system: temporary)
              bp_arr(npi)%rel_pos(1) = - body_arr(nbi)%elem(nei)%L_geom(1)/2. + body_arr(nbi)%elem(nei)%dx(1)/2. + &
                                       (i-1) * body_arr(nbi)%elem(nei)%dx(1)
              if (ncord == 2) bp_arr(npi)%rel_pos(2) = 0.
              if (ncord == 3) then
                 bp_arr(npi)%rel_pos(2) = - body_arr(nbi)%elem(nei)%L_geom(2)/2. + body_arr(nbi)%elem(nei)%dx(2)/2. + &
                              (j-1) * body_arr(nbi)%elem(nei)%dx(2)  
              endif
              bp_arr(npi)%rel_pos(3) = - body_arr(nbi)%elem(nei)%L_geom(3)/2. + body_arr(nbi)%elem(nei)%dx(3)/2. + &
                                       (k-1) * body_arr(nbi)%elem(nei)%dx(3) 
!Initial relative normal (with respect to the body axis: temporary) and area(/length in 2D) only for particles on the body surface
              bp_arr(npi)%area = 0. 
              bp_arr(npi)%normal = 0.
              if ((i==1).and.(body_arr(nbi)%elem(nei)%normal_act(5)==1)) then
                 bp_arr(npi)%normal(1) = 1.
                 if (ncord == 2) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(3)
                 if (ncord == 3) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(2)*body_arr(nbi)%elem(nei)%dx(3)
              endif
              if ((i==vec_temp(1)).and.(body_arr(nbi)%elem(nei)%normal_act(6)==1)) then
                 bp_arr(npi)%normal(1) = -1.
                 if (ncord == 2) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(3)
                 if (ncord == 3) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(2)*body_arr(nbi)%elem(nei)%dx(3)
              endif
              if (ncord == 3) then
                 if ((j==1).and.(body_arr(nbi)%elem(nei)%normal_act(3)==1)) then
                    bp_arr(npi)%normal(2) = 1.
                    bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(1)*body_arr(nbi)%elem(nei)%dx(3)
                 endif
                 if ((j==vec_temp(2)).and.(body_arr(nbi)%elem(nei)%normal_act(4)==1)) then
                    bp_arr(npi)%normal(2) = -1.
                    bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(1) * body_arr(nbi)%elem(nei)%dx(3) 
                 endif
              endif
              if ((k==1).and.(body_arr(nbi)%elem(nei)%normal_act(1)==1)) then
                 bp_arr(npi)%normal(3) = 1.
                 if (ncord == 2) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(1)
                 if (ncord == 3) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(1)*body_arr(nbi)%elem(nei)%dx(2)
              endif
              if ((k==vec_temp(3)).and.(body_arr(nbi)%elem(nei)%normal_act(2)==1)) then
                 bp_arr(npi)%normal(3) = -1.
                 if (ncord == 2) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(1)
                 if (ncord == 3) bp_arr(npi)%area = bp_arr(npi)%area + body_arr(nbi)%elem(nei)%dx(1)*body_arr(nbi)%elem(nei)%dx(2)
              endif
              mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))
              if (mod_normal > 1.) bp_arr(npi)%normal(:) = bp_arr(npi)%normal(:) / mod_normal
! Preliminary value of the relative positions (rotation of rel_pos with respect to the centre of mass of the element)
! "relative" means with respect to the element origin, in the world reference system (not in the element reference system)
              call vector_rotation(bp_arr(npi)%rel_pos,body_arr(nbi)%elem(nei)%alfa) 
              if (ncord == 2) bp_arr(npi)%rel_pos(2) = 0. 
! Preliminary value of the  absolute positions after
! translation of rel_pos of the distance corresponding element - reference system origin)
              bp_arr(npi)%pos(:) = body_arr(nbi)%elem(nei)%x_CM(:) + bp_arr(npi)%rel_pos(:)
! Body particle masses
!AA501btest start
!              xmin = - body_arr(nbi)%elem(nei)%L_geom(1)/2. + 2.*Domain%h
!              xmax = + body_arr(nbi)%elem(nei)%L_geom(1)/2. - 2.*Domain%h
!              ymin = - body_arr(nbi)%elem(nei)%L_geom(2)/2. + 2.*Domain%h
!              ymax = + body_arr(nbi)%elem(nei)%L_geom(2)/2. - 2.*Domain%h
!              zmin = - body_arr(nbi)%elem(nei)%L_geom(3)/2. + 2.*Domain%h
!              zmax = + body_arr(nbi)%elem(nei)%L_geom(3)/2. - 2.*Domain%h
!              if (ncord == 3) then
!                 if ( (bp_arr(npi)%rel_pos(1)>xmin).and.(bp_arr(npi)%rel_pos(1)<xmax) .and. &
!                      (bp_arr(npi)%rel_pos(2)>ymin).and.(bp_arr(npi)%rel_pos(2)<ymax) .and. &
!                      (bp_arr(npi)%rel_pos(3)>zmin).and.(bp_arr(npi)%rel_pos(3)<zmax) ) then 
!                     bp_arr(npi)%mass = 0.  
!                     else
!                        bp_arr(npi)%mass = body_arr(nbi)%mass / body_arr(nbi)%npart
!                 endif
!              endif
!              if (ncord == 2) then
!                 if ( (bp_arr(npi)%rel_pos(1)>xmin).and.(bp_arr(npi)%rel_pos(1)<xmax) .and. &
!                      (bp_arr(npi)%rel_pos(3)>zmin).and.(bp_arr(npi)%rel_pos(3)<zmax) ) then 
!                       bp_arr(npi)%mass = 0.  
!                     else
!                         bp_arr(npi)%mass = body_arr(nbi)%mass / body_arr(nbi)%npart
!                 endif
!              endif
!AA501btest end              
!Deactivating particle masses (before any rotation around the centre of mass of the body)
              if ( ((bp_arr(npi)%pos(1)>=body_arr(nbi)%elem(nei)%mass_deact(1)) .or. &
                   (bp_arr(npi)%pos(1)<=body_arr(nbi)%elem(nei)%mass_deact(2)) .or. &
                   (bp_arr(npi)%pos(2)>=body_arr(nbi)%elem(nei)%mass_deact(3)) .or. &
                   (bp_arr(npi)%pos(2)<=body_arr(nbi)%elem(nei)%mass_deact(4)) .or. &
                   (bp_arr(npi)%pos(3)>=body_arr(nbi)%elem(nei)%mass_deact(5)) .or. &
                   (bp_arr(npi)%pos(3)<=body_arr(nbi)%elem(nei)%mass_deact(6))) .and. &
                   (mod_normal==0.)                                                      ) then
                 bp_arr(npi)%mass = 0.
                 else
                    bp_arr(npi)%mass = body_arr(nbi)%mass / body_arr(nbi)%npart
                    aux_mass(nbi) = aux_mass(nbi) + bp_arr(npi)%mass
              endif
! Since here the relative position refers to the centre of mass of the body
              bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) - body_arr(nbi)%x_CM(:)              
! Initial normal (after rotation around the centre of mass of the element)
              call vector_rotation(bp_arr(npi)%normal,body_arr(nbi)%elem(nei)%alfa)
              if (ncord == 2) bp_arr(npi)%normal(2) = 0.
! Rotations around the centre of rotation provided in input
               bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) - body_arr(nbi)%x_rotC(:) 
               call vector_rotation(bp_arr(npi)%rel_pos,body_arr(nbi)%alfa) 
               if (ncord == 2) bp_arr(npi)%rel_pos(2) = 0.
               bp_arr(npi)%pos(:) = bp_arr(npi)%rel_pos(:) + body_arr(nbi)%x_rotC(:) 
               bp_arr(npi)%rel_pos(:) = bp_arr(npi)%pos(:) - body_arr(nbi)%x_CM(:) 
               call vector_rotation(bp_arr(npi)%normal,body_arr(nbi)%alfa)
               if (ncord == 2) bp_arr(npi)%normal(2) = 0. 
! Initial cell                                                       
              bp_arr(npi)%cell =  ParticleCellNumber(bp_arr(npi)%pos)
! Velocity
              call Vector_Product(body_arr(nbi)%omega,bp_arr(npi)%rel_pos,vec2_temp,3) 
              bp_arr(npi)%vel(:) = body_arr(nbi)%u_CM(:) + vec2_temp(:) 
              if (ncord == 2) bp_arr(npi)%vel(2) = 0. 
! Update of umax
              aux_umax = dsqrt(dot_product(bp_arr(npi)%vel,bp_arr(npi)%vel))        
              body_arr(bp_arr(npi)%body)%umax = max(body_arr(bp_arr(npi)%body)%umax,aux_umax)       

! Formal initialization of pressure and acceleration
              bp_arr(npi)%pres = 0.              
              bp_arr(npi)%acc(:) = 0.
              bp_arr(npi)%vel_mir(:) = 0.
              npi = npi + 1                                                                     
           enddo
        enddo
     enddo
!     
     enddo
  end do

!Cancel useless particles
  allocate(aux_bp_arr(n_body_part))
  do i=1,n_body_part
     if (bp_arr(i)%mass>0.) then
        aux_bp_arr(i-erased_part) = bp_arr(i)
        else
           erased_part = erased_part + 1
           body_arr(bp_arr(i)%body)%npart = body_arr(bp_arr(i)%body)%npart -1
     endif
  enddo
  n_body_part = n_body_part - erased_part 
  deallocate(bp_arr)
  allocate(bp_arr(n_body_part))
!$omp parallel do default(none) private(i) shared(n_body_part,bp_arr,aux_bp_arr)
  do i=1,n_body_part
     bp_arr(i) = aux_bp_arr(i) 
  enddo
!$omp end parallel do  

!Counting surface body particles
  do i=1,n_body_part
     if (bp_arr(i)%area>0.) n_surf_body_part = n_surf_body_part+1  
  enddo
!Allocating the list of the surface body particles
 allocate (surf_body_part(n_surf_body_part), stat = ier)
 if (ier /= 0) then
   write (nout,'(1x,a,i2)') "   surf_body_part not allocated. Error code: ",ier
 else
   write (nout,'(1x,a)') "   Arrays surf_body_part successfully allocated "
 end if
 surf_body_part = 0
!Listing surface body particles
 j = 1
 do i=1,n_body_part
    if (bp_arr(i)%area>0.) then
       surf_body_part(j) = i
       j = j+1
    endif
 enddo
 
!Redistributing lost mass 
!$omp parallel do default(none) private(i) shared(n_body_part,bp_arr,aux_mass,body_arr)
  do i=1,n_body_part
     if (bp_arr(i)%mass>0.) bp_arr(i)%mass = bp_arr(i)%mass - (aux_mass(bp_arr(i)%body)-body_arr(bp_arr(i)%body)%mass)/n_body_part
  enddo
!$omp end parallel do

! Deallocations
  deallocate(aux_bp_arr)  
!$omp parallel do default(none) private(i) shared(n_bodies,body_arr,aux_mass)
  do i=1,n_bodies
     deallocate(body_arr(i)%elem) 
  end do
!$omp end parallel do 
  deallocate(aux_mass)

  return
  end subroutine Input_Body_Dynamics
!---split


!cfile body_to_smoothing_pres.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_to_smoothing_pres
!
! Creation: Amicarelli-Agate 18October12
!
!************************************************************************************
! Module purpose : contribution of body particles to pressure smoothing
!
! Calling routines: PressureSmoothing_2D,PressureSmoothing_3D
!
! Called subroutines: /
!
! Creation: Amicarelli-Agate 13nov12
!
!************************************************************************************
!
  subroutine body_to_smoothing_pres(sompW_vec,AppUnity_vec)

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: W_vol,dis
  double precision,dimension(nag),intent(inout) :: sompW_vec,AppUnity_vec

! External functions
  double precision, external :: w

! Statements
! Loop over body particles (neighbours: fluid particles)
  do npi=1,n_body_part
! Loop over the neighbouring fluid particles 
     do j=1,nPartIntorno_bp_f(npi)
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno_bp_f(npartint)
        dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
        W_vol = w(dis,Domain%h,Domain%coefke) * ((Domain%dd/dx_dxbodies)**ncord)
        AppUnity_vec(npj) = AppUnity_vec(npj) + W_vol
        sompW_vec(npj) = sompW_vec(npj) + (bp_arr(npi)%pres-pg(npj)%pres) * W_vol
     enddo
  enddo
  
  return
  end subroutine body_to_smoothing_pres
!---split


!cfile body_to_smoothing_vel.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_to_smoothing_vel
!
!************************************************************************************
! Module purpose : contribution of body particles to pressure smoothing
!
! Creation: Amicarelli-Agate 13nov12
!
! Calling routines: inter_SmoothVelo_2D,inter_SmoothVelo_3D
!
! Called subroutines: /
!
!************************************************************************************
!
  subroutine body_to_smoothing_vel(dervel_mat,unity_vec)

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: W_vol,dis
  double precision,dimension(nag,3),intent(inout) :: dervel_mat
  double precision,dimension(nag),intent(inout) :: unity_vec

! External functions
  double precision, external :: w

! Statements
! Loop over body particles (neighbours: fluid particles)
  do npi=1,n_body_part
! Loop over the neighbouring fluid particles 
     do j=1,nPartIntorno_bp_f(npi)
        npartint = (npi-1)* NMAXPARTJ + j
        npj = PartIntorno_bp_f(npartint)
        dis = dsqrt(dot_product(rag_bp_f(:,npartint),rag_bp_f(:,npartint)))
        W_vol = w(dis,Domain%h,Domain%coefke) * ((Domain%dd/dx_dxbodies)**ncord)
        unity_vec(npj) = unity_vec(npj) + W_vol
        dervel_mat(npj,:) = dervel_mat(npj,:) + (bp_arr(npi)%vel_mir(:)-pg(npj)%vel(:)) * W_vol
     enddo
  enddo
  
  return
  end subroutine body_to_smoothing_vel
!---split



!cfile body_pressure_postpro.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: body_pressure_postpro
!
! Creation: Amicarelli-Agate 13nov12
!
!************************************************************************************
! Module purpose : Computation of pressure for body partuicles
!
! Calling routines: Loop_Irre_2D,Loop_Irre_3D
!
! Called subroutines: /
!
!************************************************************************************
!
  subroutine body_pressure_postpro

! Assigning modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE

! Declarations
  implicit none
  integer(4) :: npi,j,npartint,npj
  double precision :: Sum_W_vol,W_vol,dis,mod_normal
  double precision :: dis_vec(3)
  double precision,dimension(:),allocatable :: aux_pres
  integer,dimension(:),allocatable :: wet
  double precision :: aux

! External functions
  double precision, external :: w

! Allocations
  allocate(aux_pres(n_body_part))
  allocate(wet(n_body_part))

!Initializations
  wet = 0.

! Statements

!AA501btest start
! Loop over the body particles (neighbours: body particles at the surface of the same body)
!$omp parallel do default(none) private(Sum_W_vol,npi,npj,dis_vec,dis,mod_normal,W_vol) &
!$omp shared(aux_pres,n_body_part,bp_arr,Domain,wet)
  do npi=1,n_body_part
     aux_pres(npi) = 0.  
     mod_normal = dsqrt(dot_product(bp_arr(npi)%normal,bp_arr(npi)%normal))

!AA501btest     
!     if ((mod_normal>0.).and.(wet(npi)==1)) then
      if (mod_normal>0.) then

        Sum_W_vol = 0.
! Loop over the body particles at the surface of the same body
        do npj=1,n_body_part
           mod_normal = dsqrt(dot_product(bp_arr(npj)%normal,bp_arr(npj)%normal))
           if ( (bp_arr(npi)%body == bp_arr(npj)%body) .and. (mod_normal>0.) ) then
              dis_vec(:) = bp_arr(npj)%pos(:) - bp_arr(npi)%pos(:) 
              dis = dsqrt(dot_product(dis_vec,dis_vec))
              W_vol = w(dis,Domain%h,Domain%coefke) * bp_arr(npj)%mass
              aux_pres(npi) = aux_pres(npi) + bp_arr(npj)%pres * W_vol
              Sum_W_vol = Sum_W_vol + W_vol 
           endif
        end do
        if (Sum_W_vol > 0.000001) aux_pres(npi) = aux_pres(npi) / Sum_W_vol
     endif
  enddo
!$omp end parallel do
! Loop over the body particles: update of pressure values 
!$omp parallel do default(none) private(npi) shared(bp_arr,aux_pres,n_body_part)
  do npi=1,n_body_part

!AA501btest
     if (aux_pres(npi).ne.0.) bp_arr(npi)%pres = aux_pres(npi)

  enddo
!$omp end parallel do
!AA501btest end

! Deallocations
  deallocate(aux_pres)
  deallocate(wet)

  return
  end subroutine body_pressure_postpro
!---split

!cfile Body_dynamics_output.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Body_dynamics_output
!
!AA504 sub start
! Versions: 
! 01  Amicarelli      13nov12       (creation) .txt output for body dynamics
! 02  Amicarelli      08apr14       (v5.04) writing the surface body particle parameters      
!AA504 sub end
!
!************************************************************************************
! Module purpose : essential output for body dynamics
!
! Calling routine: Memo_Ctl
!
!************************************************************************************
!
subroutine Body_dynamics_output

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE

! Declarations ..
implicit none
integer(4)       :: nbi
character(255) :: nomefilectl_Body_dynamics
!AA504sub
integer(4)       :: npi,j,nsi
character(255) :: nomefilectl_Body_particles
double precision :: aux
! 2 auxiliary parameters: pmax_R (x>0) and pmax_L (x<0), involving boundaries pointing towards positive z
double precision, allocatable, dimension(:) :: pmax_R,pmax_L

!Allocations
allocate(pmax_R(n_bodies))
allocate(pmax_L(n_bodies))

! Loop over bodies: initialiting body maximum value of pressure (and 2 further parameters)
!$omp parallel do default(none) private(nbi) shared(n_bodies,body_arr,pmax_R,pmax_L)
 do nbi=1,n_bodies
    body_arr(nbi)%pmax = -999999.
    pmax_R(nbi) = -999999.
    pmax_L(nbi) = -999999.
 enddo
!$omp end parallel do

!Loop over body particles to estimate the 2 auxiliary parameters 
 do npi=1,n_body_part
    if (bp_arr(npi)%area > 0.) then
       body_arr(bp_arr(npi)%body)%pmax = max(body_arr(bp_arr(npi)%body)%pmax,bp_arr(npi)%pres)
       if (bp_arr(npi)%normal(3) > 0.) then
          if (bp_arr(npi)%normal(1) > 0.) then
             pmax_L(bp_arr(npi)%body) = max(pmax_L(bp_arr(npi)%body),bp_arr(npi)%pres)
             else
                pmax_R(bp_arr(npi)%body) = max(pmax_R(bp_arr(npi)%body),bp_arr(npi)%pres)
          endif     
       endif
    endif
 enddo

! File creation and heading
 write(nomefilectl_Body_dynamics,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_Body_dynamics_',it_corrente,".txt"
 open (ncpt,file=nomefilectl_Body_dynamics,status="unknown",form="formatted")
 if (it_corrente == 1) then
    write (ncpt,*) "Body dynamics values "
    write (ncpt,'(5(7x,a),3(5x,a),3(9x,a),3(1x,a),3(8x,a),(9x,a),2(7x,a))') &
        " Time(s)"," Body_ID"," x_CM(m)"," y_CM(m)"," z_CM(m)"," u_CM(m/s)"," v_CM(m/s)"," w_CM(m/s)"," Fx(N)"," Fy(N)", &
        " Fz(N)", "omega_x(rad/s)","omega_y(rad/s)","omega_z(rad/s)","Mx(N*m)","My(N*m)","Mz(N*m)","pmax(Pa)","pmax_R(Pa)", &
        "pmax_L(Pa)"
 endif
 flush(ncpt)

! Loop over the bodies
 do nbi = 1,n_bodies
    write (ncpt,'(g14.7,1x,i14,1x,18(g14.7,1x))') tempo,nbi,body_arr(nbi)%x_CM(:),body_arr(nbi)%u_CM(:),body_arr(nbi)%Force(:), &
                                               body_arr(nbi)%omega(:),body_arr(nbi)%Moment(:),body_arr(nbi)%pmax,pmax_R(nbi), &
                                               pmax_L(nbi)
 enddo
  
 close (ncpt)
 
!Monitoring the surface body particles
 write(nomefilectl_Body_particles,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_Body_particles_',it_corrente,".txt"
 open (ncpt,file=nomefilectl_Body_particles,status="unknown",form="formatted")
 if (it_corrente == 1) then
    write (ncpt,*) " Body particle parameters"
!AA601 sub
    write (ncpt,'((7x,a),(3x,a),(7x,a),3(10x,a),3(8x,a),(6x,a),4(1x,a))') &
       " Time(s)"," particle_ID"," body_ID"," x(m)"," y(m)"," z(m)"," u(m/s)"," v(m/s)"," w(m/s)"," p(N/m^2)", &
       " impact_vel_body1(m/s)"," impact_vel_body2(m/s)"," impact_vel_body3(m/s)"," impact_vel_body4(m/s)"
 endif
 flush(ncpt)
 
!AA504sub start
 do nsi=1,n_surf_body_part
    npi=surf_body_part(nsi)
    write (ncpt,'(g14.7,1x,2(i14,1x),11(g14.7,1x))') tempo,npi,bp_arr(npi)%body,bp_arr(npi)%pos(:),bp_arr(npi)%vel(:), &
    bp_arr(npi)%pres,impact_vel(nsi,1),impact_vel(nsi,2),impact_vel(nsi,3),impact_vel(nsi,4)
 enddo
!AA504sub end

deallocate(pmax_R)
deallocate(pmax_L) 
 
return
end subroutine Body_dynamics_output
!---split

!cfile Gamma_boun.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
! Creation: Amicarelli-Agate 13nov12 
!
!************************************************************************************
! Interpolative function use by Monaghan (2004, 2005) for particle-boundary interactions
double precision function Gamma_boun(r,h)
!
!.. Implicit Declarations ..
implicit none
!
!.. Local Scalars ..
double precision :: r,h,q
!
!.. Executable Statements ..
!
 q = r / h
 if ( q <= (2./3.) ) then
    Gamma_boun = 2./3.
 else if (q <= 1.0) then
    Gamma_boun = 2.*q-(3./2.)*q*q 
 else if (q <=2.0) then
    Gamma_boun = 0.5 * ((2.-q)**2)
 else
    Gamma_boun = 0.
 end if
 
 Gamma_boun = abs(Gamma_boun)

return
end function Gamma_boun 
!---split


!AA501b the whole subroutine
!AA504 sub
!cfile dis_point_plane.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : dis_point_plane
!
! Creation : Amicarelli-Agate, 25Jul12
!
!************************************************************************************
! Module purpose : Computation of the distance between a point and a plane
!AA601 sub
! Calling routine: RHS_body_dynamics,DBSPH_inlet_outlet
!
! Called routines: /
!
!************************************************************************************

subroutine dis_point_plane(P0,P1_plane,P2_plane,P3_plane,dis,normal)

! Declarations
 implicit none
 double precision,intent(IN) :: P0(3),P1_plane(3),P2_plane(3),P3_plane(3) 
 double precision,intent(INOUT) :: dis
 double precision,intent(INOUT) :: normal(3)
 double precision               :: a,b,c,d
 
! Statements
 a = (P2_plane(2)-P1_plane(2)) * (P3_plane(3)-P1_plane(3)) - (P3_plane(2)-P1_plane(2)) * (P2_plane(3)-P1_plane(3))
 b = -((P2_plane(1)-P1_plane(1)) * (P3_plane(3)-P1_plane(3)) - (P3_plane(1)-P1_plane(1)) * (P2_plane(3)-P1_plane(3)))
 c = (P2_plane(1)-P1_plane(1)) * (P3_plane(2)-P1_plane(2)) - (P3_plane(1)-P1_plane(1)) * (P2_plane(2)-P1_plane(2))
 d = -(a*P1_plane(1)+b*P1_plane(2)+c*P1_plane(3))
 normal(1) = a/dsqrt(a**2+b**2+c**2)
 normal(2) = b/dsqrt(a**2+b**2+c**2)
 normal(3) = c/dsqrt(a**2+b**2+c**2)
 dis = abs((a*P0(1)+b*P0(2)+c*P0(3)+d)/dsqrt(a**2+b**2+c**2))
 
return
end subroutine dis_point_plane
!---split



!AA501b the whole subroutine
!AA504 sub
!cfile distance_point_line_3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : distance_point_line_3D
!
! Creation : Amicarelli-Agate, 11Sep12
!
!************************************************************************************
! Module purpose : Computation of the distance between a point and a plane
!
! Calling routine: body_particles_to_continuity
!
!AA504 sub
! Called routines: three_plane_intersection
!
!************************************************************************************

subroutine distance_point_line_3D(P0,P1_line,P2_line,dis)

! Declarations
 implicit none
 double precision,intent(IN) :: P0(3),P1_line(3),P2_line(3) 
 double precision,intent(INOUT) :: dis
!AA504 
 integer(4) :: test_point
 double precision               :: a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
 double precision               :: P3_plane1(3),P3_plane2(3),intersection_pl(3),aux_vec(3)
 
! Statements 
 
!Fictitious points to define the fictitious planes whose intersection is the given line
 P3_plane1(1) = P1_line(1) - 999.
 P3_plane1(2) = P1_line(2)
 P3_plane1(3) = P1_line(3)
 P3_plane2(1) = P1_line(1) 
 P3_plane2(2) = P1_line(2) - 999.
 if (P2_line(3)==P2_line(3)) then
    P3_plane2(3) = P1_line(3) - 999.
    else
    P3_plane2(3) = P1_line(3)
 endif 
!AA504 sub comm
!Finding the coefficients for the Cartesian equation of the planes
!Plane 1: first auxiliary plane passing for the line (containing the points P1_line, P2_line and P3_plane1)
 a1 = (P2_line(2)-P1_line(2)) * (P3_plane1(3)-P1_line(3)) - (P3_plane1(2)-P1_line(2)) * (P2_line(3)-P1_line(3))
 b1 = -((P2_line(1)-P1_line(1)) * (P3_plane1(3)-P1_line(3)) - (P3_plane1(1)-P1_line(1)) * (P2_line(3)-P1_line(3)))
 c1 = (P2_line(1)-P1_line(1)) * (P3_plane1(2)-P1_line(2)) - (P3_plane1(1)-P1_line(1)) * (P2_line(2)-P1_line(2))
 d1 = -(a1*P1_line(1)+b1*P1_line(2)+c1*P1_line(3))
!AA504 sub comm
!Plane 2: second auxiliary plane passing for the line (containing the points P1_line, P2_line and P3_plane2)
 a2 = (P2_line(2)-P1_line(2)) * (P3_plane2(3)-P1_line(3)) - (P3_plane2(2)-P1_line(2)) * (P2_line(3)-P1_line(3))
 b2 = -((P2_line(1)-P1_line(1)) * (P3_plane2(3)-P1_line(3)) - (P3_plane2(1)-P1_line(1)) * (P2_line(3)-P1_line(3)))
 c2 = (P2_line(1)-P1_line(1)) * (P3_plane2(2)-P1_line(2)) - (P3_plane2(1)-P1_line(1)) * (P2_line(2)-P1_line(2))
 d2 = -(a2*P1_line(1)+b2*P1_line(2)+c2*P1_line(3))
!AA504 comm 
!Plane 3: plane passing for P0 and perpendicular to P1-P2 
 a3 = P2_line(1)-P1_line(1) 
 b3 = P2_line(2)-P1_line(2) 
 c3 = P2_line(3)-P1_line(3) 
 d3 = -(a3*P0(1)+b3*P0(2)+c3*P0(3)) 

!AA504 sub
!Intersection line-plane (intersection_pl): intersection of 3 planes
 call three_plane_intersection(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,test_point,intersection_pl)

!distance point-intersection
 aux_vec(:) = intersection_pl(:)-P0(:)
 dis = sqrt(dot_product(aux_vec,aux_vec))

 return
end subroutine distance_point_line_3D
!---split


!AA504 the whole subroutine, startig from distance_point_line_3D of AA501b
!cfile three_plane_intersection.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : three_plane_intersection
!
! Versions:
! 01   Amicarelli   08Apr14   (v5.04) New subroutine from adaptation of "distance_point_line_3D" (from v5.01, Body dynamics)
!
!************************************************************************************
! Module purpose : Computation of the intersection between 3 planes
!
! Calling routine: distance_point_line_3D,line_plane_intersection
!
! Called routines: Matrix_Inversion_3x3 
!
!************************************************************************************

subroutine three_plane_intersection(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,test_point,intersection_pl)

! Declarations
 implicit none
 double precision,intent(in)    :: a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
 double precision,intent(inout) :: intersection_pl(3)
 integer(4), intent(inout)         :: test_point
 integer(4)                     :: test
 double precision               :: b_vec(3)
 double precision               :: matrix(3,3),inverted(3,3)
 
!Statements 
!Solving the 3x3 linear system 
 matrix(1,1) = a1
 matrix(1,2) = b1
 matrix(1,3) = c1
 matrix(2,1) = a2
 matrix(2,2) = b2
 matrix(2,3) = c2
 matrix(3,1) = a3
 matrix(3,2) = b3
 matrix(3,3) = c3 
 call Matrix_Inversion_3x3(matrix,inverted,test)
 if (test==1) then
 test_point = 1
 b_vec(1) = -d1
 b_vec(2) = -d2
 b_vec(3) = -d3
 intersection_pl(1) = dot_product(inverted(1,:),b_vec) 
 intersection_pl(2) = dot_product(inverted(2,:),b_vec) 
 intersection_pl(3) = dot_product(inverted(3,:),b_vec) 
 else
 test_point = 0 
 endif
 
 return
end subroutine three_plane_intersection
!---split



!AA501b the whole subroutine
!AA504 sub
!cfile distance_point_line_2D.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : distance_point_line_2D
!
! Creation : Amicarelli, 24Oct12
!
!************************************************************************************
! Module purpose : Computation of the distance between a point and a plane
!
! Calling routine: point_inout_polygone
!
! Called routines: /
!
!************************************************************************************

subroutine distance_point_line_2D(P0,P1_line,P2_line,dis,normal)

! Declarations
 implicit none
 double precision,intent(IN) :: P0(2),P1_line(2),P2_line(2)
 double precision,intent(INOUT) :: dis
 double precision,intent(INOUT) :: normal(2)
 double precision               :: a,b,c
 
! Statements
 a = (P2_line(2)-P1_line(2))
 b = -(P2_line(1)-P1_line(1)) 
 c = -(a*P1_line(1)+b*P1_line(2))
 normal(1) = a/dsqrt(a**2+b**2)
 normal(2) = b/dsqrt(a**2+b**2)
 dis = (a*P0(1)+b*P0(2)+c)/dsqrt(a**2+b**2)
 
return
end subroutine distance_point_line_2D
!---split



!AA501b the whole subroutine
!AA504 sub
!cfile vector_rotation.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : vector_rotation
!
! Creation : Amicarelli-Agate, 2Aug12
!
!************************************************************************************
! Module purpose : Rotation of a given vector provided the vector of the rotation 
!                  angles (3D)
!
! Calling routine: time_integration_body_dynamics,Input_Body_Dynamics
!
! Called routines: MatrixProduct
!
!************************************************************************************

subroutine vector_rotation(vector,angle)

! Declarations
 implicit none
 double precision,intent(IN) :: angle(3) 
 double precision,intent(INOUT) :: vector(3)
 double precision               :: vec_temp(3)
 double precision               :: mat1_temp(3,3),mat2_temp(3,3),mat3_temp(3,3),mat4_temp(3,3),cos_dir(3,3)
 
! Statements
  mat1_temp(1,1) = 1.
  mat1_temp(1,2) = 0.
  mat1_temp(1,3) = 0.
  mat1_temp(2,1) = 0.
  mat1_temp(2,2) = cos(angle(1))
  mat1_temp(2,3) = sin(-angle(1))
  mat1_temp(3,1) = 0.
  mat1_temp(3,2) = -sin(-angle(1))
  mat1_temp(3,3) = cos(-angle(1))
  mat2_temp(1,1) = cos(-angle(2))
  mat2_temp(1,2) = 0.
  mat2_temp(1,3) = -sin(-angle(2))
  mat2_temp(2,1) = 0. 
  mat2_temp(2,2) = 1.
  mat2_temp(2,3) = 0.
  mat2_temp(3,1) = sin(-angle(2)) 
  mat2_temp(3,2) = 0.
  mat2_temp(3,3) = cos(-angle(2))
  mat3_temp(1,1) = cos(-angle(3))
  mat3_temp(1,2) = sin(-angle(3))
  mat3_temp(1,3) = 0.
  mat3_temp(2,1) = -sin(-angle(3))
  mat3_temp(2,2) = cos(-angle(3))
  mat3_temp(2,3) = 0.
  mat3_temp(3,1) = 0. 
  mat3_temp(3,2) = 0.
  mat3_temp(3,3) = 1.     
  call MatrixProduct(mat1_temp,mat2_temp,mat4_temp,3,3,3)
  call MatrixProduct(mat4_temp,mat3_temp,cos_dir,3,3,3)
  vec_temp(:) = vector(:) 
  vector(1) = dot_product(cos_dir(1,:),vec_temp)
  vector(2) = dot_product(cos_dir(2,:),vec_temp)
  vector(3) = dot_product(cos_dir(3,:),vec_temp)

return
end subroutine vector_rotation
!---split



!AA501b the whole subroutine
!AA504 sub
!cfile reference_system_change.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : reference_system_change
!
! Creation : Amicarelli-Agate, 24Oct12
!
!************************************************************************************
! Module purpose : Transformation of the coordinates, expressed in a new reference system
!
!AA601
! Calling routine: RHS_body_dynamics,DBSPH_inlet_outlet
!
! Called routines: MatrixTransposition 
!
!************************************************************************************

subroutine reference_system_change(pos_old_ref,new_origin_old_ref,new_cos_dir_old_ref,pos_new_ref)

! Declarations
 implicit none
 integer(4) :: i,aux_int
 double precision,intent(IN) :: pos_old_ref(3),new_origin_old_ref(3) 
 double precision,intent(IN) :: new_cos_dir_old_ref(3,3)
 double precision,intent(INOUT) :: pos_new_ref(3)
 double precision :: aux_vec(3)
 double precision :: inv_new_cos_dir_old_ref(3,3)
 
! Statements
 aux_vec(:) = pos_old_ref(:) - new_origin_old_ref(:)
 call Matrix_Inversion_3x3(new_cos_dir_old_ref,inv_new_cos_dir_old_ref,aux_int)
 do i=1,3
    pos_new_ref(i) = dot_product(inv_new_cos_dir_old_ref(i,:),aux_vec)
 enddo

return
end subroutine reference_system_change
!---split



!AA501b the whole subroutine
!AA504 sub
!cfile Body_dynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : point_inout_polygone
!
!AA504 sub start
! Versions: 
! 01  Amicarelli      24Oct12       (creation) 
! 02  Amicarelli      08apr14       (v5.04) removed minor errors      
!AA504 sub end
!
!************************************************************************************
! Module purpose : test to evaluate if a point lies inside or strictly outside a polygone
!                  (triangle or quadrilateral)
!
!AA504sub
!AA601sub
! Calling routine: RHS_body_dynamics,GeneratePart,DBSPH_inlet_outlet
!
! Called routines: distance_point_line_2D
!
!************************************************************************************

subroutine point_inout_polygone(point,n_sides,point_pol_1,point_pol_2,point_pol_3,point_pol_4,test)

! Declarations
 implicit none
 integer(4),intent(in) :: n_sides
 double precision,intent(in) :: point(2),point_pol_1(2),point_pol_2(2),point_pol_3(2),point_pol_4(2)
 integer(4),intent(inout) :: test
 double precision :: dis1,dis2
 double precision :: normal(2)
 
! Statements
  call distance_point_line_2D(point,point_pol_1,point_pol_2,dis1,normal)
!AA504 rm line 
  if (dis1.ne.0.) call distance_point_line_2D(point_pol_3,point_pol_1,point_pol_2,dis2,normal)
  if (dsign(dis1,dis2).ne.(dis1)) then
     test = 0
     else
        call distance_point_line_2D(point,point_pol_2,point_pol_3,dis1,normal)
!AA504 rm line
        if (dis1.ne.0.) call distance_point_line_2D(point_pol_1,point_pol_2,point_pol_3,dis2,normal)
        if (dsign(dis1,dis2).ne.(dis1)) then
           test = 0
           else
              if (n_sides == 3) then
                 call distance_point_line_2D(point,point_pol_3,point_pol_1,dis1,normal)
                 if (dis1.ne.0.) call distance_point_line_2D(point_pol_2,point_pol_3,point_pol_1,dis2,normal) 
                 if (dsign(dis1,dis2).ne.(dis1)) then
                    test = 0
                    else
                    test = 1
                 endif
                 else
                    call distance_point_line_2D(point,point_pol_3,point_pol_4,dis1,normal)
                    if (dis1.ne.0.) call distance_point_line_2D(point_pol_1,point_pol_3,point_pol_4,dis2,normal) 
                    if (dsign(dis1,dis2).ne.(dis1)) then
                       test = 0
                       else
                          call distance_point_line_2D(point,point_pol_4,point_pol_1,dis1,normal)
                          if (dis1.ne.0.) call distance_point_line_2D(point_pol_2,point_pol_4,point_pol_1,dis2,normal) 
                          if (dsign(dis1,dis2).ne.(dis1)) then
                             test = 0
                             else
                                test = 1
                          endif
                    endif
              endif
        endif
  endif
        
return
end subroutine point_inout_polygone
!---split

!AA501b the whole subroutine
!AA504 sub
!cfile Body_dynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : Matrix_Inversion_3x3
!
!AA504 sub start
! Versions: 
! 01  Amicarelli      19Jul12       (creation) 
! 02  Amicarelli      08apr14       (v5.04) treatment of bad conditioned matrices      
!AA504 sub end
!
!************************************************************************************
! Module purpose : computation of the inverse (inv) of a provided 3x3 matrix (mat)
!
! Calling routine: RHS_body_dynamics
!
! Called routines: /
!
!************************************************************************************

!AA504 sub
subroutine Matrix_Inversion_3x3(mat,inv,test)

! Declarations
 implicit none
 double precision,intent(IN) :: mat(3,3) 
 double precision,intent(INOUT) :: inv(3,3)
!AA504
! integer(4),intent(inout),optional :: test  
 integer(4),intent(inout) :: test  
 double precision :: det

! Statements
 det = mat(1,1)*mat(2,2)*mat(3,3) + mat(2,1)*mat(3,2)*mat(1,3) + mat(3,1)*mat(1,2)*mat(2,3) &
      -mat(1,1)*mat(3,2)*mat(2,3) - mat(3,1)*mat(2,2)*mat(1,3) - mat(2,1)*mat(1,2)*mat(3,3)
!AA504 start
 if (det.ne.0) then
 test = 1
!AA504 end
 inv(1,1) = mat(2,2)*mat(3,3) - mat(2,3)*mat(3,2)
 inv(1,2) = mat(1,3)*mat(3,2) - mat(1,2)*mat(3,3)
 inv(1,3) = mat(1,2)*mat(2,3) - mat(1,3)*mat(2,2)
 inv(2,1) = mat(2,3)*mat(3,1) - mat(2,1)*mat(3,3)
 inv(2,2) = mat(1,1)*mat(3,3) - mat(1,3)*mat(3,1)
 inv(2,3) = mat(1,3)*mat(2,1) - mat(1,1)*mat(2,3)
 inv(3,1) = mat(2,1)*mat(3,2) - mat(2,2)*mat(3,1)
 inv(3,2) = mat(1,2)*mat(3,1) - mat(1,1)*mat(3,2)
 inv(3,3) = mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)
 inv(:,:) = inv(:,:) / det
!AA504
 else
 test = 0
 endif

return
end subroutine Matrix_Inversion_3x3
!---split

!AA501b the whole subroutine
!AA504 sub
!cfile Body_dynamics.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : Matrix_Inversion_2x2
!
! Creation : Amicarelli-Agate, 19 July 2012
!
!************************************************************************************
! Module purpose : computation of the inverse (inv) of a provided 3x3 matrix (mat)
!
! Calling routine: RHS_body_dynamics
!
! Called routines: /
!
!************************************************************************************

subroutine Matrix_Inversion_2x2(mat,inv)

! Declarations
 implicit none
 double precision,intent(IN) :: mat(2,2) 
 double precision,intent(INOUT) :: inv(2,2)
 double precision :: det

! Statements
 det = mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2)
 inv(1,1) = mat(2,2)
 inv(1,2) = -mat(2,1)
 inv(2,1) = -mat(1,2)
 inv(2,2) = mat(1,1)
 inv(:,:) = inv(:,:) / det

return
end subroutine Matrix_Inversion_2x2
!---split
