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

