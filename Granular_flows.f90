!cfile mixture_viscosity.f90
!AA504 all the subroutines

!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: mixture_viscosity
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose :  (v5.04) To compute the mixture viscosity (granular flows)
!
! Calling routines: Loop_Irre_3D,Loop_Irre_2D
!
! Called subroutines: / 
!
!************************************************************************************

  subroutine mixture_viscosity

! Assigning modules
  use GLOBAL_MODULE 
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use FILES_ENTITIES

! Declarations
  implicit none
  double precision :: lambda_star,z_int,mu_main_fluid,vel2,tau_c,mu_e,vel,lambda 
  integer(4) :: npi,i_cell,i_aux,i_grid,j_grid,k_grid
  integer(4),external :: ParticleCellNumber,CellIndices 
  
! Statements
  lambda_star = 0.0001d0
  mu_main_fluid = Med(Granular_flows_options%ID_main_fluid)%visc * Med(Granular_flows_options%ID_main_fluid)%den0
!$omp parallel do default(none) shared(nag,Med,pg,lambda_star,Granular_flows_options,mu_main_fluid,Domain,ind_interfaces,nout,BoundaryDataPointer) &
!$omp private(npi,i_cell,i_aux,i_grid,j_grid,k_grid,z_int,vel2,lambda,mu_e,tau_c,vel)
 do npi=1,nag
     pg(npi)%Bn = 0.d0
     if ((Med(pg(npi)%imed)%tipo=="granular").and.(pg(npi)%state=="flu")) then
        i_cell = ParticleCellNumber(pg(npi)%coord) 
        i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
!AA601 sub
        if ( (Granular_flows_options%viscosity_blt_formula==1) .or. (Granular_flows_options%viscosity_blt_formula==4) ) then
           if (ind_interfaces(i_grid,j_grid,3).ne.0) then                  
              z_int = pg(ind_interfaces(i_grid,j_grid,3))%coord(3)
              pg(npi)%sigma_prime = (-Domain%grav(3)) * &
                                    (Med(Granular_flows_options%ID_granular)%den0-Med(Granular_flows_options%ID_main_fluid)%den0) * &
                                    (z_int-pg(npi)%coord(3))
              else
! this case (probably) never occurs
                 pg(npi)%sigma_prime = pg(npi)%pres * (Med(Granular_flows_options%ID_granular)%den0-Med(Granular_flows_options%ID_main_fluid)%den0) &
                                       / Med(Granular_flows_options%ID_granular)%den0             
           endif 
           if (pg(npi)%sigma_prime<0.d0) pg(npi)%sigma_prime = 0.d0
           if (Granular_flows_options%viscosity_blt_formula==3) pg(npi)%sigma_prime = 0.d0
!    secinv=sqrt(I2(e_ij) (secinv is the sqrt of the second inveriant of the strain rate tensor)
           tau_c = pg(npi)%sigma_prime*dtan(Med(Granular_flows_options%ID_granular)%phi)
           mu_e = mu_main_fluid * (1.0d0+5.0d0/2.0d0*Med(pg(npi)%imed)%gran_vol_frac_max)
!AA601 sub           
           if (Granular_flows_options%viscosity_blt_formula==1) then 
              lambda = lambda_star
              else           
                 vel = dsqrt(dot_product(pg(npi)%vel,pg(npi)%vel))
                 if (vel>0.d0) then
                    pg(npi)%Bn = tau_c*(z_int-pg(npi)%coord(3))/(vel*mu_e)
                    else
                       pg(npi)%Bn = 1000.d0
                 endif
                 if (pg(npi)%Bn>1000.d0) pg(npi)%Bn = 1000.d0    
                 lambda = min(lambda_star,lambda_star/pg(npi)%Bn)
           endif      
           pg(npi)%mu = mu_e + tau_c / (lambda+pg(npi)%secinv)
           elseif (Granular_flows_options%viscosity_blt_formula==2) then
              lambda = lambda_star  
              vel2 = dot_product(pg(npi)%vel,pg(npi)%vel)
              pg(npi)%mu = mu_main_fluid + Granular_flows_options%Chezy_friction_coeff * &
                           Med(Granular_flows_options%ID_main_fluid)%den0 * vel2 / (lambda+pg(npi)%secinv) 
              elseif (Granular_flows_options%viscosity_blt_formula==3) then   
                 mu_e = mu_main_fluid * (1.0d0+5.0d0/2.0d0*Med(pg(npi)%imed)%gran_vol_frac_max)
                 pg(npi)%mu = mu_e
        endif    
        pg(npi)%visc = pg(npi)%mu/pg(npi)%dens
     endif
 end do
!$omp end parallel do

 return
 end subroutine mixture_viscosity
!---split



!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: compute_k_BetaGamma
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) To compute k_BetaGamma=teta_c/teta_c,0
!                  k_BetaGamma is the ratio between Shields critical non-dimensional 
!                  stress for a generic 3D slope (teta_c) and its analogous value 
!                  defined by Shields' diagram (teta_c,00) on a flat bed
!
! Calling routines: Shields
!
! Called subroutines: quadratic_equation 
!
!************************************************************************************

  subroutine compute_k_BetaGamma(npi,npj,DistZmin)

! Assigning modules
  use GLOBAL_MODULE 
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use files_entities

! Declarations
  implicit none
  integer(4), intent(in) :: npi,npj
  double precision,intent(in) :: DistZmin
   integer(4) :: i_cell,i_aux,i_grid,j_grid,k_grid,n_roots,test_root1,test_root2,case_test,inpj_grid,jnpj_grid
  double precision :: k_v,k_Beta0,U_inf,Re,DELTA_Seminara,aux_var,a,b,c,root1,root2,Uf,z_aux,alfa
  double precision :: U_inf_vec(2),Uf_vec(2)
  double precision :: s_uni_vec(3),n2_uni_vec(3),vel_n(3)
  integer(4), external :: CellIndices,ParticleCellNumber
  
! Initializations
  k_v = 0.41d0  
  test_root1 = 0
  test_root2 = 0
  
! Statements
! unity vector aligned with the component of the trajectory normal to the interface
  vel_n(:) = dot_product(pg(npj)%vel_old,pg(npi)%normal_int) * pg(npi)%normal_int(:)  
  s_uni_vec(:) = pg(npj)%vel_old(:) - vel_n(:)
  s_uni_vec(:) = s_uni_vec(:)/dsqrt(dot_product(s_uni_vec,s_uni_vec)) 
! n2 is the bi-normal to the particle trajectory     
  call Vector_Product(pg(npi)%normal_int,s_uni_vec,n2_uni_vec,3)
! Beta
  if(ncord==3) then
     pg(npi)%Beta_slope = dasin(-s_uni_vec(3))
     else
        alfa = dasin(pg(npi)%normal_int(3))
        pg(npi)%Beta_slope = PIGRECO/2.d0 - alfa
        if (s_uni_vec(3)>0.d0) pg(npi)%Beta_slope = -pg(npi)%Beta_slope
  endif
! Gamma     
  pg(npi)%Gamma_slope = dasin(dabs(n2_uni_vec(3)))
  if (Granular_flows_options%Gamma_slope_flag==0) pg(npi)%Gamma_slope = 0.d0 
  if ( (pg(npi)%Beta_slope>=Med(pg(npi)%imed)%Phi) .or. (dabs(pg(npi)%Gamma_slope)>=Med(pg(npi)%imed)%Phi) ) then
! If the internal friction angle is equal or less than a slope angle (Beta or Gamma) then the granular particle is eroded
        pg(npi)%k_BetaGamma = 0.0d0
        else
! Compute k_Beta0
           k_Beta0 = dcos(pg(npi)%Beta_slope) * (1-dtan(pg(npi)%Beta_slope)/dtan(Med(pg(npi)%imed)%Phi))
           if (pg(npi)%Gamma_slope==0.d0) then
! compute k_BetaGamma=k_Beta0
              pg(npi)%k_BetaGamma = k_Beta0
              else
! Getting U_inf
                 i_cell = ParticleCellNumber(pg(npi)%coord) 
                 i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
                 U_inf = dsqrt(dot_product(pg(npj)%vel_old,pg(npj)%vel_old)) 
! Compute Re (flow around a sphere)
                 Re = U_inf*Med(pg(npi)%imed)%D50/pg(npj)%visc 
! Compute the drag coefficient C_D (Morrison 2013)              
                 if (Re<=100.0d0) then
                    pg(npi)%C_D = 1.0d0
                    else
                       if (Re>=1.0e6) then
                          pg(npi)%C_D = 0.2d0
                          else
                             pg(npi)%C_D = 24.0d0/Re + (2.6d0*Re/5.0d0)/(1.0d0+(Re/5.0d0)**1.52d0) +  &
                                   (0.411d0*(Re/2.63e5)**(-7.94d0))/(1.0d0+(Re/2.63e5)**(-8.0d0)) + (Re**0.8d0)/4.61e5
                       endif
                 endif    
! Compute the lift coefficient C_L (our regression from data in Seminara et al. 2002 WRR)                
                 if (Re<=2000.) then
                    pg(npi)%C_L = 0.07d0
                    else
                       if (Re>=17.5e3) then
                          pg(npi)%C_L = 0.5d0
                          else
                             pg(npi)%C_L = (9.0d0*(10.0d0)**(-5.0d0))*(Re**0.882d0)
                       endif         
                 endif    
! Compute DELTA_Seminara
                 Uf = U_inf
! Compute height over fluid - bed load interface
                 z_aux = DistZmin
                 DELTA_Seminara = (4.0d0/3.0d0) * dtan(Med(pg(npi)%imed)%Phi) * (pg(npi)%C_L/pg(npi)%C_D) * (pg(npi)%u_star*Med(pg(npi)%imed)%D50)/(k_v*Uf*z_aux)   
                 if (pg(npi)%Beta_slope==0.0d0) then
! Compute k_BetaGamma=k_0Gamma
                      pg(npi)%k_BetaGamma = dcos(pg(npi)%Gamma_slope) * (dsqrt(1.0d0 - ((1.0d0-DELTA_Seminara)**2.0d0) * &
                                       (dtan(pg(npi)%Gamma_slope)/dtan(Med(pg(npi)%imed)%Phi))**2.0d0) - DELTA_Seminara) / (1-DELTA_Seminara)
                    else
! compute coefficients (a,b,c) of the quadratic equation for k_BetaGamma and solve it
                       aux_var = (dtan(pg(npi)%Beta_slope))**2.0d0 + (dtan(pg(npi)%Gamma_slope))**2.0d0
                       a = (1.0d0-DELTA_Seminara)
                       b = 2.0d0 * ( DELTA_Seminara/dsqrt(1.0d0+aux_var) + dsin(pg(npi)%Beta_slope)/dtan(Med(pg(npi)%imed)%Phi) )
                       c = (1.0d0+DELTA_Seminara) / (1.0d0+aux_var) * (-1.0d0+aux_var/(dtan(Med(pg(npi)%imed)%Phi))**2.0d0)                    
                       call quadratic_equation(a,b,c,n_roots,root1,root2)
                       if ((root1<=k_Beta0).and.(root1>=0.d0)) test_root1 = 1
                       if ((root2<=k_Beta0).and.(root2>=0.d0)) test_root2 = 1
                       case_test = n_roots * 100 + test_root1*10 + test_root2
                       select case (case_test)
                          case(211)
                                pg(npi)%k_BetaGamma = max(root1,root2)
                          case(210,110)
                             pg(npi)%k_BetaGamma = root1 
                          case(201)
                             pg(npi)%k_BetaGamma = root2
                          case(200,100,0)
                                pg(npi)%k_BetaGamma = k_Beta0   
                       end select
                    endif
          endif
    endif
  
  return
  end subroutine compute_k_BetaGamma
!---split



!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: quadratic_equation
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) To solve a quadratic equation compute k_s 
!
! Calling routines: compute_k_BetaGamma
!
! Called subroutines: / 
!
!************************************************************************************

 subroutine quadratic_equation(a,b,c,n_roots,root1,root2)

! Assigning modules
 use GLOBAL_MODULE 
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations
 implicit none
 double precision, intent(in) :: a,b,c
 integer(4),intent(out) :: n_roots
 double precision, intent(out) :: root1,root2
 double precision :: discriminant
 
! Statements
 root1=-999.d0
 root2=-999.d0
 if (a.ne.0.d0) then
    discriminant = b**2 - 4.0d0*a*c
    if (discriminant>0.0d0) then
       n_roots = 2
       root1 = (-b-dsqrt(discriminant))/(2.0d0*a)
       root2 = (-b+dsqrt(discriminant))/(2.0d0*a)
       else
       if (discriminant==0.0d0) then
          n_roots = 1
          root1 = -b/(2.0d0*a)
          else
             n_roots = 0
       endif      
    endif
    else
       if (b.ne.0.d0) then
          n_roots = 1
          root1 = -c/b
       endif   
  endif       
  
 return
 end subroutine quadratic_equation
!---split



!cfile write_Granular_flows_interfaces.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : write_Granular_flows_interfaces 
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Writiing of the interfaces (bed load transport / granular flows)
!
! Calling routine: Loop_Irre_3D/Loop_Irre_2D 
!
! Called subroutine: /
!
!************************************************************************************

subroutine write_Granular_flows_interfaces

! Modules
 use FILES_ENTITIES
 use GLOBAL_MODULE
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations ..
 implicit none
 integer(4)     :: i_grid,j_grid,i_cell,i_aux,i,k_grid,i2_grid,j2_grid
 character(255) :: nomefile_blt_interfaces
 double precision :: x_grid,y_grid,z_free_surface,z_BedLoad_PureFluid,z_bed
 double precision :: pos_aux(3)
 integer(4), external :: CellIndices,ParticleCellNumber

! “.txt” file creation and heading
 write(nomefile_blt_interfaces,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_blt_interfaces_',it_corrente,".txt"
 open (ncpt,file=nomefile_blt_interfaces,status="unknown",form="formatted")
 
 if (it_corrente == 1) then
! First step
    write (ncpt,*) "Bed load transport interfaces "
    write (ncpt,'((7x,a),(5x,a),(5x,a),(7x,a),(7x,a),(8x,a))') &
           " Time(s)"," x_grid(m)"," y_grid(m)"," z_FS(m)"," z_fm(m)"," z_b(m)"
    flush(ncpt)
    else
! Other steps 
! Loop over the monitoring lines
       do i=1,Granular_flows_options%monitoring_lines
          if (Granular_flows_options%lines(i,1)==-999.d0) then
              do i2_grid=1,Grid%ncd(1)
                 pos_aux(1) = (i2_grid-0.5) * Grid%dcd(1) + grid%extr(1,1)
                 pos_aux(2) = Granular_flows_options%lines(i,2)
                 pos_aux(3) = 0.5 * Grid%dcd(3) + grid%extr(3,1) !pos_aux(3) = 0.d0
                 i_cell = ParticleCellNumber(pos_aux) 
                 i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
                 if (ind_interfaces(i_grid,j_grid,1)>0) then
                    z_free_surface = pg(ind_interfaces(i_grid,j_grid,1))%coord(3) + Domain%dd/2.d0
                    else
                         z_free_surface = -999.d0
                 endif
                 if (ind_interfaces(i_grid,j_grid,4)>0) then
                    z_bed = pg(ind_interfaces(i_grid,j_grid,4))%coord(3) + Domain%dd/2.d0
                    else
                         z_bed = -999.d0
                 endif 
                    if (ind_interfaces(i_grid,j_grid,3)>0) then
                    z_BedLoad_PureFluid = pg(ind_interfaces(i_grid,j_grid,3))%coord(3) + Domain%dd/2.d0
                    else
                         z_BedLoad_PureFluid = z_bed 
                 endif   
                 write(ncpt,'(6(g14.6,1x))') tempo,pos_aux(1),pos_aux(2),z_free_surface,z_BedLoad_PureFluid,z_bed
              end do 
          endif 
          if (Granular_flows_options%lines(i,2)==-999.d0) then
              do j2_grid=1,Grid%ncd(2)
                 pos_aux(1) = Granular_flows_options%lines(i,1)
                 pos_aux(2) = (j2_grid-0.5) * Grid%dcd(2) + grid%extr(2,1)
                 pos_aux(3) = 0.5 * Grid%dcd(3) + grid%extr(3,1) !pos_aux(3) = 0.d0
                 i_cell = ParticleCellNumber(pos_aux) 
                 i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
                 if (ind_interfaces(i_grid,j_grid,1)>0) then
                    z_free_surface = pg(ind_interfaces(i_grid,j_grid,1))%coord(3) + Domain%dd/2.d0
                    else
                         z_free_surface = -999.d0
                 endif
                 if (ind_interfaces(i_grid,j_grid,4)>0) then
                    z_bed = pg(ind_interfaces(i_grid,j_grid,4))%coord(3) + Domain%dd/2.d0
                    else
                         z_bed = -999.d0
                 endif 
                    if (ind_interfaces(i_grid,j_grid,3)>0) then
                    z_BedLoad_PureFluid = pg(ind_interfaces(i_grid,j_grid,3))%coord(3) + Domain%dd/2.d0
                    else
                         z_BedLoad_PureFluid = z_bed 
                 endif   
                 write(ncpt,'(6(g14.6,1x))') tempo,pos_aux(1),pos_aux(2),z_free_surface,z_BedLoad_PureFluid,z_bed
              end do
          endif
          if ((Granular_flows_options%lines(i,1).ne.-999.d0).and.(Granular_flows_options%lines(i,2).ne.-999.d0)) then
              pos_aux(1) = Granular_flows_options%lines(i,1)
              pos_aux(2) = Granular_flows_options%lines(i,2)
              pos_aux(3) = 0.5 * Grid%dcd(3) + grid%extr(3,1) !pos_aux(3) = 0.d0
              i_cell = ParticleCellNumber(pos_aux) 
              i_aux = CellIndices(i_cell,i_grid,j_grid,k_grid)
              if (ind_interfaces(i_grid,j_grid,1)>0) then
                 z_free_surface = pg(ind_interfaces(i_grid,j_grid,1))%coord(3) + Domain%dd/2.d0
                 else
                      z_free_surface = -999.d0
              endif
              if (ind_interfaces(i_grid,j_grid,4)>0) then
                 z_bed = pg(ind_interfaces(i_grid,j_grid,4))%coord(3) + Domain%dd/2.d0
                 else
                      z_bed = -999.d0
              endif 
                 if (ind_interfaces(i_grid,j_grid,3)>0) then
                 z_BedLoad_PureFluid = pg(ind_interfaces(i_grid,j_grid,3))%coord(3) + Domain%dd/2.d0
                 else
                      z_BedLoad_PureFluid = z_bed 
              endif   
              write(ncpt,'(6(g14.6,1x))') tempo,pos_aux(1),pos_aux(2),z_free_surface,z_BedLoad_PureFluid,z_bed
          endif          
       end do
 endif 
 
 close (ncpt)
 
 return
end subroutine write_Granular_flows_interfaces
!---split



!cfile Update_Zmax_at_grid_vert_columns.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : Update_Zmax_at_grid_vert_columns
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Updating the 2D array of the maximum values of the fluid particle
!                  height, at the nodes of the grid vertical columns (only in 3D).
!                  Printing the 2D field of the water depth (current time step), 
!                  according to output frequency chosen in input file (only in 3D).
!                  Printing the 2D fields of specific flow rate components (current
!                  time step), at the same frequency of the water depth (only in 3D).
!
! Calling routine: Loop_Irre_3D
!
! Called routines: /
!
!************************************************************************************

subroutine Update_Zmax_at_grid_vert_columns(print_flag)

!Using modules
 use GLOBAL_MODULE
 use AdM_USER_TYPE
 use ALLOC_MODULE
 use files_entities

! Declarations
 implicit none
 integer(4) :: npi,GridColumn
 double precision :: pos(3)
 double precision, allocatable, dimension(:) :: Z_fluid_step,h_step,qx_step,qy_step,qx_step_grid,qy_step_grid,n_part_step 
 integer(4) :: i_zone,i_vertex,i_aux,i_grid,j_grid
 integer(4), intent(in) :: print_flag
 character(255) :: nomefile_h_step
 
!External function 
 integer(4),external :: ParticleCellNumber

!Allocations
 do i_zone=1,NPartZone
    if (Partz(i_zone)%IC_source_type == 2) then
       allocate(h_step(Partz(i_zone)%npoints))
       allocate(qx_step(Partz(i_zone)%npoints))
       allocate(qy_step(Partz(i_zone)%npoints))
       exit
    endif
 enddo
 allocate(Z_fluid_step(Grid%ncd(1)*Grid%ncd(2)))
 allocate(qx_step_grid(Grid%ncd(1)*Grid%ncd(2)))
 allocate(qy_step_grid(Grid%ncd(1)*Grid%ncd(2)))
 allocate(n_part_step(Grid%ncd(1)*Grid%ncd(2)))
 Z_fluid_step = -999.d0
 h_step = 0.d0
 qx_step = 0.d0
 qy_step = 0.d0
 qx_step_grid = 0.d0
 qy_step_grid = 0.d0 
 n_part_step = 0
 
! Statements
!$omp parallel do default(none) shared(nag,pg,Grid,Z_fluid_max,Domain,Z_fluid_step,n_part_step,qx_step_grid,qy_step_grid,tempo) private(npi,GridColumn,pos)
 do npi=1,nag
    pos(1) = pg(npi)%coord(1)
    pos(2) = pg(npi)%coord(2)
    pos(3) = Grid%extr(3,1) + 0.0000001
    GridColumn = ParticleCellNumber(pos)
    Z_fluid_step(GridColumn) = max(Z_fluid_step(GridColumn),(pg(npi)%coord(3)+0.5*Domain%dd)) 
    qx_step_grid(GridColumn) = qx_step_grid(GridColumn) + pg(npi)%vel(1)
    qy_step_grid(GridColumn) = qy_step_grid(GridColumn) + pg(npi)%vel(2)
    n_part_step(GridColumn) = n_part_step(GridColumn) + 1
    Z_fluid_max(GridColumn) = max (Z_fluid_max(GridColumn),(pg(npi)%coord(3)+0.5*Domain%dd))
 end do
!$omp end parallel do
!$omp parallel do default(none) shared(Grid,n_part_step,qx_step_grid,qy_step_grid,tempo) private(GridColumn,pos,i_grid,j_grid)
    do i_grid=1,Grid%ncd(1)
       do j_grid=1,Grid%ncd(2)
          pos(1) = (i_grid-0.5) * Grid%dcd(1) + Grid%extr(1,1)
          pos(2) = (j_grid-0.5) * Grid%dcd(2) + Grid%extr(2,1)
          pos(3) = Grid%extr(3,1) + 0.0000001
          GridColumn = ParticleCellNumber(pos)
          if (n_part_step(GridColumn)>0) then
             qx_step_grid(GridColumn) = qx_step_grid(GridColumn)/n_part_step(GridColumn)
             qy_step_grid(GridColumn) = qy_step_grid(GridColumn)/n_part_step(GridColumn)
             else
                qx_step_grid(GridColumn) = 0.d0
                qy_step_grid(GridColumn) = 0.d0
          endif
       enddo
    enddo
!$omp end parallel do
! “.txt” file creation and heading (only at the beginning of the simulation)
 if (it_corrente == 1) then
    write(nomefile_h_step,"(a,a)") nomecaso(1:len_trim(nomecaso)),"_h_qx_qy_step.txt"
    open(ncpt,file=nomefile_h_step,status="unknown",form="formatted")
    write(ncpt,*) "Water depth (h), Specific flow rate components (q_x,q_y)"
    write(ncpt,'(7(a))') "           x(m)","           y(m)","  h_step(m)"," Z_fl_max_stp(m)","     z_topog(m)","     q_x(m^2/s)","     q_y(m^2/s)"
    flush(ncpt)
    else 
! Writing the 2D free surface field at the current time step
       if (print_flag==1) then
          write(nomefile_h_step,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_h_qx_qy_step',it_corrente,".txt"
          open (ncpt,file=nomefile_h_step,status="unknown",form="formatted")
       endif   
       do i_zone=1,NPartZone
          if (Partz(i_zone)%IC_source_type == 2) then
!$omp parallel do default(none) shared(ncpt,Partz,Vertice,Grid,h_step,Z_fluid_step,i_zone,qx_step,qy_step,qx_step_grid,qy_step_grid,n_part_step,q_max,tempo,print_flag) &
!$omp private(i_vertex,GridColumn,pos,i_aux)
             do i_vertex=Partz(i_zone)%ID_first_vertex,Partz(i_zone)%ID_last_vertex
                i_aux = i_vertex-Partz(i_zone)%ID_first_vertex+1 
                pos(1) =  Vertice(1,i_vertex)
                pos(2) = Vertice(2,i_vertex)
                pos(3) = Grid%extr(3,1)+0.0000001
                GridColumn = ParticleCellNumber(pos)
                h_step(i_aux) = max ((Z_fluid_step(GridColumn)-Vertice(3,i_vertex)),0.d0)
                qx_step(i_aux) = qx_step_grid(GridColumn)
                qy_step(i_aux) = qy_step_grid(GridColumn)
                if (qx_step_grid(GridColumn)/=-999.d0) then
                    qx_step(i_aux) = qx_step_grid(GridColumn)*h_step(i_aux)
                    qy_step(i_aux) = qy_step_grid(GridColumn)*h_step(i_aux)
                    q_max(i_aux) = max(q_max(i_aux),dsqrt(qx_step(i_aux)**2 + qy_step(i_aux)**2))
                endif    
                if (print_flag==1) then
!$omp critical (omp_write_h_step)
                   write (ncpt,'(7(f14.4,1x))') Vertice(1,i_vertex),Vertice(2,i_vertex),h_step(i_aux),Z_fluid_step(GridColumn),Vertice(3,i_vertex), &
                                                qx_step(i_aux),qy_step(i_aux)
!$omp end critical (omp_write_h_step)
                endif
             end do
!$omp end parallel do    
             exit
          endif
       enddo
 endif

! Closing the step file
 if (print_flag==1) close (ncpt)
 
!Deallocations
 if (allocated(Z_fluid_step)) deallocate(Z_fluid_step)
 if (allocated(h_step)) deallocate(h_step)
 if (allocated(qx_step)) deallocate(qx_step)
 if (allocated(qy_step)) deallocate(qy_step)
 if (allocated(qx_step_grid)) deallocate(qx_step_grid)
 if (allocated(qy_step_grid)) deallocate(qy_step_grid)
 if (allocated(n_part_step)) deallocate(n_part_step)
 
 return
end subroutine Update_Zmax_at_grid_vert_columns
!---split



!cfile write_h_max.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : write_h_max
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Compute and write the 2D array of the maximum values of the water
!                  depth, at the nodes of a Cartesian topography (only in 3D),
!                  together with the 2D field of the maximum (over time) specific 
!                  flow rates.
!
! Calling routine: Gest_Trans
!
! Called routines: /
!
!************************************************************************************

subroutine write_h_max

!Using modules
 use GLOBAL_MODULE
 use AdM_USER_TYPE
 use ALLOC_MODULE
 use files_entities
 
!Declarations
 implicit none
 integer(4) :: i_zone,i_vertex,GridColumn
 double precision :: pos(3)
 double precision,dimension(:),allocatable :: h_max 
 character(255) :: nomefile_h_max
 
!External function
 integer(4),external :: ParticleCellNumber
 
! h_max .txt file: creation and heading 
 write(nomefile_h_max,"(a,a)") nomecaso(1:len_trim(nomecaso)),"_h_max_q_max.txt"
 open(ncpt,file=nomefile_h_max,status="unknown",form="formatted")
 write(ncpt,*) "Maximum water depth (m) and specific flow rate (m^2/s)"
 write(ncpt,'(6(a))') "           x(m)","           y(m)","       h_max(m)"," Z_fluid_max(m)","     z_topog(m)","   q_max(m^2/s)"
 flush(ncpt) 
 
!Statements
 do i_zone=1,NPartZone
    if (Partz(i_zone)%IC_source_type == 2) then
!Allocating h_max 
       allocate(h_max(Partz(i_zone)%npoints))
!Initializing h_max
       h_max = 0.
!$omp parallel do default(none) shared(Partz,Vertice,Grid,h_max,Z_fluid_max,ncpt,i_zone,q_max) private(i_vertex,GridColumn,pos)
       do i_vertex=Partz(i_zone)%ID_first_vertex,Partz(i_zone)%ID_last_vertex
          pos(1) = Vertice(1,i_vertex)
          pos(2) = Vertice(2,i_vertex)
          pos(3) = Grid%extr(3,1)+0.0000001
          GridColumn = ParticleCellNumber(pos)
          h_max(i_vertex-Partz(i_zone)%ID_first_vertex+1) = max ((Z_fluid_max(GridColumn)-Vertice(3,i_vertex)),0.d0)
!$omp critical (omp_write_h_max)
          write(ncpt,'(6(f14.4,1x))')Vertice(1,i_vertex),Vertice(2,i_vertex),h_max(i_vertex-Partz(i_zone)%ID_first_vertex+1),Z_fluid_max(GridColumn),Vertice(3,i_vertex) &
                                     ,q_max(i_vertex-Partz(i_zone)%ID_first_vertex+1)        
!$omp end critical (omp_write_h_max)
       end do
!$omp end parallel do    
       exit
    endif
 end do

! h_max .txt file: closing
 close (ncpt)

!Deallocations
 if (allocated(h_max)) deallocate(h_max)
 if (allocated(q_max)) deallocate(q_max)

 return
end subroutine write_h_max
!---split



!cfile sub_Q_sections.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : sub_Q_sections 
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Writing flow rate at several sections (only in 3D)
!
! Calling routine: Loop_Irre_3D 
!
! Called subroutine: line_plane_intersection,point_inout_polygone,Vector_Product,
!                    reference_system_change
!
!************************************************************************************

subroutine sub_Q_sections

! Modules
 use FILES_ENTITIES
 use GLOBAL_MODULE
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations ..
 implicit none
 integer(4)     :: npi,i_sect,test_intersection_point,test_inout,j
 character(255) :: nomefile_Q_sections
 double precision :: sign_dis_part_sect_old,sign_dis_part_sect
 double precision :: P_intersection(3),aux_vec_1(3),aux_vec_2(3),P_intersection_loc(3)
 integer(4), allocatable, dimension(:) :: n_particles 

!Allocations
 allocate(n_particles(Q_sections%n_sect))
 
! Initializing n_particles
 n_particles = 0

! “.txt” file creation and heading (only at the beginning of the simulation)
 write(nomefile_Q_sections,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_Q_sections_',it_corrente,".txt"
 open (ncpt,file=nomefile_Q_sections,status="unknown",form="formatted")
 
 if (it_corrente == 1) then
! First step
    write (ncpt,*) "Flow rate values(m^3/s) "
    write (ncpt,'((7x,a),(4x,a),(5x,a),(6x,a),(5x,a),(3x,a),(5x,a),6(5x,a))') &
           " Time(s)"," ID_section","fluid_type"," Q(m^3/s)"," Area(m^2)"," n_particles"," dt_out(s)"," norm_x(m)"," norm_y(m)"," norm_z(m)"," x1_sec(m)"," y1_sec(m)"," z1_sec(m)"
    flush(ncpt)
!Initializing fluid particle old coordinates 
!$omp parallel do default(none) shared(nag,Q_sections,pg) private(npi)
 do npi=1,nag
    pg(npi)%sect_old_pos(:) = pg(npi)%coord(:)
 end do  
!$omp end parallel do
!Initializing section local axis and vertex local coordinates
!$omp parallel do default(none) shared(Q_sections) private(i_sect)
 do i_sect=1,Q_sections%n_sect
!Local axis 1: Vertex2-vertex1
    Q_sections%section(i_sect)%loc_axis(:,1) = Q_sections%section(i_sect)%vertex(2,:)-Q_sections%section(i_sect)%vertex(1,:)
    Q_sections%section(i_sect)%loc_axis(:,1) = Q_sections%section(i_sect)%loc_axis(:,1) / &
                                               dsqrt(dot_product(Q_sections%section(i_sect)%loc_axis(:,1),Q_sections%section(i_sect)%loc_axis(:,1)))
!Local axis 2: opposite of the vector product of Local axis 1 and normal 
    call Vector_Product(Q_sections%section(i_sect)%loc_axis(:,1),Q_sections%section(i_sect)%normal,Q_sections%section(i_sect)%loc_axis(:,2),3)
    Q_sections%section(i_sect)%loc_axis(:,2) = -Q_sections%section(i_sect)%loc_axis(:,2)
!Local axis 3: normal
    Q_sections%section(i_sect)%loc_axis(:,3) = Q_sections%section(i_sect)%normal(:)
    Q_sections%section(i_sect)%vertex_loc(1,:) = 0.0d0
    call reference_system_change(Q_sections%section(i_sect)%vertex(2,:),Q_sections%section(i_sect)%vertex(1,:), &
                                 Q_sections%section(i_sect)%loc_axis,Q_sections%section(i_sect)%vertex_loc(2,:)) 
    call reference_system_change(Q_sections%section(i_sect)%vertex(3,:),Q_sections%section(i_sect)%vertex(1,:), &
                                 Q_sections%section(i_sect)%loc_axis,Q_sections%section(i_sect)%vertex_loc(3,:)) 
    if (Q_sections%section(i_sect)%n_vertices==4) then
        call reference_system_change(Q_sections%section(i_sect)%vertex(4,:),Q_sections%section(i_sect)%vertex(1,:), &
                                     Q_sections%section(i_sect)%loc_axis,Q_sections%section(i_sect)%vertex_loc(4,:))
    endif
! Allocating and initializing the flow rate vector
    allocate(Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1))
    Q_sections%section(i_sect)%flow_rate(:) = 0.0d0
 end do
 !$omp end parallel do 
 else
!Other steps     
!Loop over the fluid particles
!$omp parallel do default(none) shared(nag,Q_sections,pg,n_particles,Med) &
!$omp private(npi,i_sect,sign_dis_part_sect_old,sign_dis_part_sect,test_intersection_point,P_intersection,test_inout,aux_vec_1,aux_vec_2,P_intersection_loc)
 do npi=1,nag
!Loop over the flow rate monitoring sections 
    do i_sect=1,Q_sections%n_sect
!Check if a fluid particle has crossed the plane of the section
       aux_vec_1 = pg(npi)%sect_old_pos - Q_sections%section(i_sect)%vertex(1,:)
       aux_vec_2 = pg(npi)%coord - Q_sections%section(i_sect)%vertex(1,:)
       sign_dis_part_sect_old = dot_product(aux_vec_1,Q_sections%section(i_sect)%normal)
       sign_dis_part_sect = dot_product(aux_vec_2,Q_sections%section(i_sect)%normal)
       if ((sign_dis_part_sect*sign_dis_part_sect_old)<0) then
!Check if the intersection between the line (defined by the old and new particle position) and the plane (containing the section) is a point            
          call line_plane_intersection(pg(npi)%sect_old_pos,pg(npi)%coord,Q_sections%section(i_sect)%vertex(1,:),Q_sections%section(i_sect)%vertex(2,:), &
                                       Q_sections%section(i_sect)%vertex(3,:),test_intersection_point,P_intersection)
          if (test_intersection_point==1) then
!Local coordinates
             call reference_system_change(P_intersection,Q_sections%section(i_sect)%vertex(1,:),Q_sections%section(i_sect)%loc_axis,P_intersection_loc) 
!Check if the intersection point lies within the section
             call point_inout_polygone(P_intersection_loc,Q_sections%section(i_sect)%n_vertices,Q_sections%section(i_sect)%vertex_loc(1,:),Q_sections%section(i_sect)%vertex_loc(2,:), &
                                       Q_sections%section(i_sect)%vertex_loc(3,:),Q_sections%section(i_sect)%vertex_loc(4,:),test_inout)
             if (test_inout==1) then
!$omp critical (omp_Q_Sections)                 
!Update the number of crossed particles and the flowing mass                  
                n_particles(i_sect) = n_particles(i_sect) + 1
                Q_sections%section(i_sect)%flow_rate(pg(npi)%imed) = Q_sections%section(i_sect)%flow_rate(pg(npi)%imed) + pg(npi)%mass/Med(pg(npi)%imed)%den0 &
                                                                     *  (sign_dis_part_sect)/abs(sign_dis_part_sect)
!$omp end critical (omp_Q_Sections)
             endif
          endif
       endif
    end do
! update old coordinates    
    pg(npi)%sect_old_pos(:) = pg(npi)%coord(:)
 end do  
!$omp end parallel do 
!Loop over the flow rate monitoring sections 
 do i_sect=1,Q_sections%n_sect
    do j=1,Q_sections%n_fluid_types
!Compute the volume to flow rate per fluid type: flowing mass by time 
       Q_sections%section(i_sect)%flow_rate(j) = Q_sections%section(i_sect)%flow_rate(j)/Q_sections%dt_out
!Writing the flow rate on a ".txt" file    
       write(ncpt,'((f14.6,1x),2(i14,1x),2(f14.6,1x),(i14,1x),7(f14.6,1x))') &
             tempo,i_sect,j,Q_sections%section(i_sect)%flow_rate(j),Q_sections%section(i_sect)%area,n_particles(i_sect),Q_sections%dt_out,Q_sections%section(i_sect)%normal(:), &
             Q_sections%section(i_sect)%vertex(1,:)
!Update the global volume flow rate: flowing mass by time 
       Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1) = Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1) + Q_sections%section(i_sect)%flow_rate(j)       
!Zeroing the flow rate per fluid type    
       Q_sections%section(i_sect)%flow_rate(j) = 0.0d0
    end do
!Writing the global flow rate on a ".txt" file  
    j = Q_sections%n_fluid_types + 1
    write(ncpt,'((f14.6,1x),2(i14,1x),2(f14.6,1x),(i14,1x),7(f14.6,1x))') &
          tempo,i_sect,j,Q_sections%section(i_sect)%flow_rate(j),Q_sections%section(i_sect)%area,n_particles(i_sect),Q_sections%dt_out,Q_sections%section(i_sect)%normal(:), &
          Q_sections%section(i_sect)%vertex(1,:) 
!Zeroing the flow rate per fluid type    
    Q_sections%section(i_sect)%flow_rate(Q_sections%n_fluid_types+1) = 0.0d0    
 end do
 
 endif 
!endif on first step 
 
 close (ncpt)
 
!Deallocations
 deallocate(n_particles) 
 
 return
end subroutine sub_Q_sections
!---split



!cfile line_plane_intersection.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : line_plane_intersection
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation, adapting the subroutine distance_point_line_3D (from v5.02, Body_dynamics)
!
!************************************************************************************
! Module purpose : (v5.04) Computation of the intersection point, if unique, between a line and a plane
!
! Calling routine: sub_Q_sections
!
! Called routines: three_plane_intersection
!
!************************************************************************************

subroutine line_plane_intersection(P1_line,P2_line,P1_plane3,P2_plane3,P3_plane3,test_point,intersection_pl)

! Declarations
 implicit none
 double precision,intent(IN) :: P1_line(3),P2_line(3),P1_plane3(3),P2_plane3(3),P3_plane3(3) 
 integer(4),intent(inout) :: test_point
 double precision,intent(INOUT) :: intersection_pl(3)
 double precision               :: a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3
 double precision               :: P3_plane1(3),P3_plane2(3),b_vec(3)
 
!Statements 
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
!Finding the coefficients for the Cartesian equation of the planes
!Plane 1: first auxiliary plane passing for the line (containing the points P1_line, P2_line and P3_plane1)
 a1 = (P2_line(2)-P1_line(2)) * (P3_plane1(3)-P1_line(3)) - (P3_plane1(2)-P1_line(2)) * (P2_line(3)-P1_line(3))
 b1 = -((P2_line(1)-P1_line(1)) * (P3_plane1(3)-P1_line(3)) - (P3_plane1(1)-P1_line(1)) * (P2_line(3)-P1_line(3)))
 c1 = (P2_line(1)-P1_line(1)) * (P3_plane1(2)-P1_line(2)) - (P3_plane1(1)-P1_line(1)) * (P2_line(2)-P1_line(2))
 d1 = -(a1*P1_line(1)+b1*P1_line(2)+c1*P1_line(3))
!Plane 2: second auxiliary plane passing for the line (containing the points P1_line, P2_line and P3_plane2)
 a2 = (P2_line(2)-P1_line(2)) * (P3_plane2(3)-P1_line(3)) - (P3_plane2(2)-P1_line(2)) * (P2_line(3)-P1_line(3))
 b2 = -((P2_line(1)-P1_line(1)) * (P3_plane2(3)-P1_line(3)) - (P3_plane2(1)-P1_line(1)) * (P2_line(3)-P1_line(3)))
 c2 = (P2_line(1)-P1_line(1)) * (P3_plane2(2)-P1_line(2)) - (P3_plane2(1)-P1_line(1)) * (P2_line(2)-P1_line(2))
 d2 = -(a2*P1_line(1)+b2*P1_line(2)+c2*P1_line(3))
!Plane 3: plane passing for P1_plane3, P2_plane3 and P3_plane3
 a3 = (P2_plane3(2)-P1_plane3(2)) * (P3_plane3(3)-P1_plane3(3)) - (P3_plane3(2)-P1_plane3(2)) * (P2_plane3(3)-P1_plane3(3))
 b3 = -((P2_plane3(1)-P1_plane3(1)) * (P3_plane3(3)-P1_plane3(3)) - (P3_plane3(1)-P1_plane3(1)) * (P2_plane3(3)-P1_plane3(3)))
 c3 = (P2_plane3(1)-P1_plane3(1)) * (P3_plane3(2)-P1_plane3(2)) - (P3_plane3(1)-P1_plane3(1)) * (P2_plane3(2)-P1_plane3(2))
 d3 = -(a3*P1_plane3(1)+b3*P1_plane3(2)+c3*P1_plane3(3))

!Intersection line-plane (intersection_pl): intersection of 3 planes
 call three_plane_intersection(a1,b1,c1,d1,a2,b2,c2,d2,a3,b3,c3,d3,test_point,intersection_pl)
 
 return
end subroutine line_plane_intersection
!---split



!AA601 the whole subroutine
!cfile Granular_flows.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : area_quadrilateral
!
! Versions: 
! 01 Amicarelli 26Jan15 (creation, AA601, DBSPH-input)
!
!************************************************************************************
! Module purpose : Computation of the area of a generic quadrilateral from the coordinates of its vertices
!
! Calling routine: ReadSectionFlowRate
!
! Called routines: area_triangle
!
!************************************************************************************

subroutine area_quadrilateral(P1,P2,P3,P4,area)

! Declarations
 implicit none
 double precision,intent(IN)    :: P1(3),P2(3),P3(3),P4(3)
 double precision,intent(INOUT) :: area
 double precision               :: area_triangle_1,area_triangle_2
!AA601
 double precision               :: aux_normal(3)
  
!Statements 
!Area of the triangle P1,P2,P3: 0.5*vector_product(vec_a1,vec_b), vec_a1=(P2-P1),vec_b=(P3-P1)
 call area_triangle(P1,P2,P3,area_triangle_1,aux_normal)
!Area of the trinagle P1,P4,P3: 0.5*vector_product(vec_a2,vec_b), vec_a2=(P4-P1),vec_b=(P3-P1)
 call area_triangle(P1,P4,P3,area_triangle_2,aux_normal)
!Area of the quadrilateral: sum of the areas of the trinagles
 area = area_triangle_1 + area_triangle_2
 return
end subroutine area_quadrilateral
!---split



!AA601 the whole subroutine
!cfile Granular_flows.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : area_triangle
!
! Creation        : Amicarelli A., 31July14
!
!************************************************************************************
! Module purpose : Computation of the area of a generic triangle, 
!                  provided the coordinates of its vertices
!
! Calling routine: area_quadrilateral,Import_ply_surface_meshes
!
! Called routines: Vector_Product
!
!************************************************************************************

subroutine area_triangle(P1,P2,P3,area,normal)

! Declarations
 implicit none
 double precision,intent(IN)    :: P1(3),P2(3),P3(3)
 double precision,intent(OUT)   :: area
 double precision,intent(OUT)   :: normal(3)
 double precision               :: vec_1(3),vec_2(3)

 !Statements 
 vec_1(:) = P2(:)-P1(:)
 vec_2(:) = P3(:)-P1(:) 
 call Vector_Product(vec_1,vec_2,normal,3)
 area = 0.5d0 * dsqrt(dot_product(normal,normal))
 normal(:) = normal(:)/(2.d0*area)

 return
end subroutine area_triangle
!---split



!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: initialization_fixed_granular_particle
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) To initialize the most of the fixed SPH granular particles 
!
! Calling routines: Shields,SetParticleParameters,Loop_Irre_2D,Loop_Irre_3D
!
! Called subroutines: / 
!
!************************************************************************************

 subroutine initialization_fixed_granular_particle(npi)

! Assigning modules
 use GLOBAL_MODULE 
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations
 implicit none
 integer(4),intent(in)  :: npi
 
! Statements
 pg(npi)%Beta_slope = -999.d0
 pg(npi)%Gamma_slope = -999.d0 
 pg(npi)%u_star = 0.d0
 pg(npi)%C_L = 0.d0
 pg(npi)%C_D = 0.d0
 pg(npi)%k_BetaGamma = -999.d0
 pg(npi)%tau_tauc = 0.0d0
 
 return
 end subroutine initialization_fixed_granular_particle
!---split



!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name: fixed_bed_slope_limited
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Forced deposition or no erosion for particles at least 2h below 
!                  the fixed (as it is defined in the column9 during the same time step: i.e. the maximum 
!                  slope of the fixed bed is 2h/2h. This avoid an eventual too fast propagation of erosion within the column (erosion is a surface phenomena).
!
! Calling routines: Shields
!
! Called subroutines: / 
!
!************************************************************************************

 subroutine fixed_bed_slope_limited(npi,igridi,jgridi,test)

! Assigning modules
 use GLOBAL_MODULE 
 use AdM_USER_TYPE
 use ALLOC_MODULE

! Declarations
 implicit none
 integer(4),intent(in) :: npi,igridi,jgridi
 logical,intent(inout) :: test
 double precision :: Velocity2,fixed_bed_tolerance
 double precision :: aux_vec(3)
 integer(4) :: aux_ID
 double precision :: pretot
!Initializations 
 if (ncord==3) then
     fixed_bed_tolerance = 4.d0*Domain%h
     elseif (ncord==2) then
        fixed_bed_tolerance = 2.d0*Domain%h
 endif    
! Statements
 if (ind_interfaces(igridi,jgridi,4)>0) then
    if (pg(npi)%coord(3)<(pg(ind_interfaces(igridi,jgridi,4))%coord(3)-fixed_bed_tolerance)) then
       if (pg(npi)%state=="flu") then
          pg(npi)%state = "sol"   
          pg(npi)%vel = 0.d0
          pg(npi)%var = 0.d0
          pg(npi)%sigma_prime = 0.0d0  
       endif
       if (pg(npi)%indneighliqsol.ne.0) then
          aux_ID = pg(npi)%indneighliqsol
          elseif (pg(npi)%state=="flu") then 
             aux_ID = pg(npi)%ind_neigh_mob_for_granmob    
             else
                aux_ID = pg(npi)%ind_neigh_mix_bed 
       endif
       if (aux_ID.ne.0) then
          aux_vec(:) = pg(aux_ID)%vel_old(:) - pg(npi)%vel_old(:)
          Velocity2 = dot_product(aux_vec,aux_vec) 
          pretot = pg(aux_ID)%pres  + (pg(aux_ID)%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(pg(aux_ID)%imed)%den0 + Velocity2 * half * Med(pg(aux_ID)%imed)%den0
          pg(npi)%dens = med(pg(npi)%imed)%den0 + (pretot / (Med(pg(npi)%imed)%celerita*Med(pg(npi)%imed)%celerita))   
          else
             if (ind_interfaces(igridi,jgridi,3)>0) then
                pretot = pg(ind_interfaces(igridi,jgridi,3))%pres  + &
                         (pg(ind_interfaces(igridi,jgridi,3))%coord(3) - pg(npi)%coord(3)) * (-Domain%grav(3)) * med(pg(npi)%imed)%den0
                pg(npi)%dens = med(pg(npi)%imed)%den0 + (pretot / (Med(pg(npi)%imed)%celerita*Med(pg(npi)%imed)%celerita))
             endif
       endif 
    test = .false.
    endif 
 endif

 return
 end subroutine fixed_bed_slope_limited
!---split


