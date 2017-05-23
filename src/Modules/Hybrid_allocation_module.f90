!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: Hybrid_allocation_module            
! Description: Module to define derived types of both dynamically and statically
!              allocated variables (Di Monaco et al., 2011, EACFM; Manenti et 
!              al., 2012; JHE; Amicarelli et al., 2013, IJNME; Amicarelli et 
!              al., 2015, CAF).                    
!-------------------------------------------------------------------------------
module Hybrid_allocation_module
use Static_allocation_module
type TyGlobal
   logical          :: NormFix                        
   logical          :: Slip
   integer(4)       :: iplot_fr ! Output frequency on log file                      
   integer(4)       :: imemo_fr ! Saving frequency                      
   integer(4)       :: irest_fr ! Restart file frequency                      
   integer(4)       :: icpoi_fr ! Monitoring point output frequency 
   integer(4)       :: ipllb_fr ! Output frequency for free surface                      
   integer(4)       :: ipllb_md ! Reference fluid for free surface
   integer(4)       :: istart ! Initial time step ID                         
   integer(4)       :: ioutopt ! Printing code for log file                       
   integer(4)       :: itmax ! Maximum step number   
   integer(4)       :: body_part_reorder ! Flag for body particles reordering
   integer(4)       :: nag_aux ! Rough and slight overestimation of the number 
                               ! of fluid particles in the reservoir (auxiliary 
                               ! parameter useful for extruded reservoir IC)           
   integer(4)       :: MAXCLOSEBOUNDFACES ! Maximum number of neighbouring 
                                          ! SA-SPH faces for a computational 
                                          ! particle
   integer(4)       :: MAXNUMCONVEXEDGES ! Maximum number of convex edges             
   integer(4)       :: density_thresholds ! Density_thresholds flag (default=0;
                                          ! =1 for very low bulk modulus  
                                          ! -preliminary simulations-)            
   integer(4)       :: time_split ! Flag for Leapfrog scheme                    
   integer(4)       :: RKscheme ! RK scheme (1,2,3,4) (time_split=0)                      
   integer(4)       :: time_stage ! Stage of RK schemes                                           
   double precision :: tmax ! Maximum physical time                           
   double precision :: dx ! dx: particle size                             
   double precision :: trunc ! h/dx, h=dx*trunc                         
   double precision :: PVolume ! dx^D, D: domain dimensionality                        
   double precision :: coefke  
   double precision :: coefkacl  
   double precision :: CFL ! CFL number
   double precision :: vsc_coeff ! viscous stability condition coefficient, 
                                 ! default value: 0.05
   double precision :: prif ! Reference pressure                          
   double precision :: plot_fr ! Output frequency for log file                       
   double precision :: memo_fr ! Frequency for result saving 
   double precision :: rest_fr ! Output frequency for restart files
   double precision :: cpoi_fr ! Output frequency for monitoring points
   double precision :: pllb_fr ! Output frequency for free surface                       
   double precision :: depth_dt_out ! Output frequency for water depth                 
   double precision :: depth_it_out_last ! Last time step when water depth 
                                         ! was printed (auxiliary variable)             
   double precision :: TetaP ! Partial smoothing parameter for pressure                         
   double precision :: TetaV ! Partial smoothing parameter for velocity                         
   double precision :: h ! Kernel support length scale                              
   double precision :: start ! Simulation start time                         
   double precision :: COEFNMAXPARTI ! Max number of fluid 
                                     ! particles=nag*COEFNMAXPARTI                 
   double precision :: COEFNMAXPARTJ ! maxb(maximum number of 
                                     ! fluid neighbours) = COEFNMAXPARTJ * 
                                     ! (4h/dx)^D  
   double precision :: t0 ! Time at the beginning of the simulation 
                          ! (origin: beginning of the year)                            
   double precision :: t_pre_iter ! Time at the beginning of the iterations
                                  ! (origin: beginning of the year)                 
   double precision :: grav(3) ! Gravity                        
   double precision :: coord(3,2) ! Coordinates of 2 vertices of a diagonal 
                                  ! of the parallelepiped domain                       
   character(4)     :: tipo
   character(100)    :: file
   character(1)     :: Psurf 
   character(1)     :: RandomPos ! IC particle distribution noise. "r": slight 
                                 ! white noise is added, otherwise nothing.  
end type TyGlobal

! Background positioning grid 
type TyGriglia
   integer(4)       :: nmax ! Number of cells                           
   integer(4)       :: ncd(3) ! Number of cells for each axis direction                         
   double precision :: dcd(3) ! Grid resolutions                          
   double precision :: extr(3,2) ! Coordinates of 2 vertices of a diagonal 
                                 ! of the parallelepiped domain                      
end type TyGriglia

! Fluid particle 
type TyParticle
   integer(4)       :: DBSPH_inlet_ID ! ID of an inlet neighbouring DB-SPH 
                                      ! surface element, if any (otherwise it is
                                      ! null)                
   integer(4)       :: DBSPH_outlet_ID ! ID of an outlet neighbouring DB-SPH  
                                       ! surface element, if any (otherwise it
                                       ! is null)               
   integer(4)       :: indneighliqsol ! Neighbouring "liquid"/"granular" 
                                      ! particle for
                                      ! a computational "granular"/"liquid" 
                                      ! particle (erosion criterion)                
   integer(4)       :: ind_neigh_mix_bed ! Index of the mobile ("flu"/mixture) 
                                         ! particle, which is the closest to a 
                                         ! computational fixed ("sol"/bed) 
                                         ! particle and viceversa             
   integer(4)       :: blt_flag ! Bed-load transport flag: 0:generic position,
                                ! 1:free surface; 2:top of bed-load transport
                                ! zone, 3:fixed bed                      
   integer(4)       :: ind_neigh_mob_for_granmob ! Index of the mobile 
                                                 ! ("flu"/mixture) particle,
                                                 ! which is the closest to
                                                 ! the computational mobile
                                                 ! granular ("sol"/bed) one     
   integer(4)       :: kodvel ! Velocity code                         
   integer(4)       :: koddens ! Density code                       
   integer(4)       :: CloseBcOut                     
   integer(4)       :: cella ! Cell number                          
   integer(4)       :: izona                         
   integer(4)       :: icol ! Colour 
   integer(4)       :: imed ! Fluid ID                          
   integer(4)       :: FS ! Free Surface (0,3: no free surface; 1,2: free 
                          ! surface) (2,3: deactivated for Shep evolution) 
   integer(4)       :: laminar_flag ! Laminar flag
   double precision :: densass                        
   double precision :: mass ! Mass                          
   double precision :: dens ! Density                          
   double precision :: pres ! Pressure                           
   double precision :: dden ! Continuity equation LHS                             
   double precision :: uni ! Shepard coefficient (Shep: discrete or integral,
                           ! it depends on the input data)                           
   double precision :: vpres ! Pressure correction                          
   double precision :: vden ! Density correction                           
   double precision :: secinv ! Second invariant of the strain rate tensor 
                              ! for incompressible fluids                        
   double precision :: dudx  
   double precision :: dudy 
   double precision :: dvdx 
   double precision :: dvdy 
   double precision :: visc ! Kinematic viscosity                          
   double precision :: mu ! Dynamic viscosity; equivalent mixture dynamic
                          ! viscosity in case of bed-laod transport layer                             
   double precision :: tstop ! Stop time                          
   double precision :: mno                       
   double precision :: ang                       
   double precision :: VolFra                         
   double precision :: rhoc                           
   double precision :: rhow                          
   double precision :: tiroc                          
   double precision :: cden                           
   double precision :: wden                           
   double precision :: diffu                         
   double precision :: coefdif ! Diffusion coefficient                        
   double precision :: IntEn ! Specific internal energy                         
   double precision :: Envar ! Partial smoothing contributions
                             ! for specific internal energy                          
   double precision :: dEdT ! Energy equation LHS                          
   double precision :: Csound ! Sound speed                        
   double precision :: DensShep ! Density times Shepard coefficient                      
   double precision :: rhoSPH_new ! SPH approximation of density at the  
                                  ! on-going time step 
   double precision :: rhoSPH_old ! SPH approximation of density at the 
                                  ! previous time step 
   double precision :: dShep ! Lagrangian derivative of Shepard 
                             ! coefficient 
   double precision :: sigma ! Discrete Shepard coefficient
   double precision :: sigma_same_fluid ! Discrete Shepard coefficient involving
                                        ! neighbours of the same fluid   
   double precision :: Gamma ! Integral Shepard coefficient
   double precision :: Gamma_last_active ! Last value of Gamma before FS=3
   double precision :: dens_init_err ! Initial difference between SPH 
                                     ! approx. of density and its exact 
                                     ! value. This is added in the interior
                                     ! domain to each SPH density approx. 
                                     ! to solve inlet problems.
   double precision :: Beta_slope ! Main slope angle of the fixed bed 
                                  ! (along the direction aligned with 
                                  ! the mean flow; bed-load transport)
   double precision :: Gamma_slope ! Transversal slope angle of the fixed bed
                                   ! (along the direction transversal to the 
                                   ! mean flow; bed-load transport)
   double precision :: sigma_prime_m ! Mean of the effective normal stresses
   double precision :: pres_fluid ! Pressure of the fluid phase (bed-load 
                                  ! transport)
   double precision :: u_star ! Friction velocity representative of the
                              ! Surface Neutral Boundary Layer 
   double precision :: C_L ! Lift coefficient for 3D erosion criterion
   double precision :: C_D ! Drag coefficient for 3D erosion criterion
   double precision :: tau_tauc ! Ratio (tau/tau_c) between bottom shear
                                ! stress and critical shear stress (from
                                ! erosion criterion)
   double precision :: k_BetaGamma ! Ratio between 3D critical shear stress
                                   ! (Shields-Seminara) and analogous 2D 
                                   ! criterion (tauc/tauc,00) 
   double precision :: normal_int(3) ! Normal of the interface between the
                                     ! mobile particles and the fixed particles 
                                     ! (bed-load transport in the presence of an
                                     ! erosion criterion); the normal vector 
                                     ! points intward the mobile domain
   double precision :: normal_int_mixture_top(3) ! Normal of the interface 
                                                 ! between the mixture and the 
                                                 ! pure fluid (top of the 
                                                 ! bed-load transport layer, 
                                                 ! mixture side)
   double precision :: normal_int_sat_top(3) ! Normal of the interface between 
                                             ! the fully saturated mixture and 
                                             ! the rest of the domain
   double precision :: vel_old(3) ! Velocity vector before the erosion 
                                  ! criterion (3D erosion criterion)
   double precision :: normal_int_old(3) ! normal_int variable before the 
                                         ! erosion criterion (3D erosion 
                                         ! criterion)
   double precision :: drho(3) ! Density gradient (SPH pseudo-consistent 
                               ! approximation over fluid particles)
   double precision :: veldif(3) ! Velocity in diffusion model 
   double precision :: rijtempmin(3) ! Minimum distance between a SPH mixture
                                     ! particle and a SPH pure fluid particle
   double precision :: vstart(3) ! Initial velocity for fixed particles 
   double precision :: velmorr(3) ! Velocity to use for Morris scheme
   double precision :: zer(3) 
   double precision :: var(3) ! Partially smoothed velocity 
   double precision :: acc(3) ! Acceleration
   double precision :: coord(3) ! Position
   double precision :: CoordOld(3) ! Position at the previous time step 
   double precision :: sect_old_pos(3) ! Position at the
                                       ! previous flow rate writing step   
   double precision :: vel(3) ! Velocity
   double precision :: velass(3) ! Imposed velocity 
   double precision :: dvel(3,3) ! Velocity gradient (SPH pseudo-consistent
                                 ! approximation over fluid particles)
   character(1)     :: slip ! BC conditions for fixed particles (f=free slip;
                            ! n=no slip; c=cont slip)

   character(3)     :: vel_type ! Movement type
   character(3)     :: state ! Particle "status" ("sol": fixed; "flu": mobile)
end type TyParticle

! Particle intermediate (time stage) values (time integration scheme RK2)
type Tytime_stage
   double precision :: ts_dden ! Stage continuity equation LHS 
   double precision :: ts_dens ! Stage density
   double precision :: ts_IntEn ! Stage specific internal energy 
   double precision :: ts_dEdT ! Stage energy equation LHS
   double precision :: ts_coord(3) ! Stage position 
   double precision :: ts_var(3) ! Stage partially smoothed velocity
   double precision :: ts_vel(3) ! Stage velocity
   double precision :: ts_acc(3) ! Stage acceleration
end type Tytime_stage

! Semi-particles and wall elements (DB-SPH)
type TyParticle_w
   integer(4)       :: cella ! Cell number 
   integer(4)       :: adjacent_faces(3) ! Index of the adjacent_faces (3 
                                         ! at maximum)
   integer(4)       :: wet ! =1 for wet wall particle (distance from fluid 
                           ! particle smaller than 1.3dx)
   integer(4)       :: surface_mesh_file_ID ! ID of the surface mesh file
   double precision :: dens ! Density
   double precision :: pres ! Pressure
   double precision :: weight ! Area(3D)/length(2D) of the wall element 
   double precision :: mass ! Mass of the semi-particle
   double precision :: k_d ! Depth coefficient
   double precision :: volume ! Semi-particle volume (area in 2D)
   double precision :: sigma ! Discrete Shepard coefficient of the wall elements
                             ! depending on fluid particles (not on
                             ! semi-particles)  
   double precision :: kin_visc_semi_part ! Kinematic viscosity of the
                                          ! semi-particle                             
   double precision :: normal(3) ! Normal
   double precision :: coord(3) ! Position
   double precision :: vel(3) ! Velocity
   double precision :: grad_vel_VSL_times_mu(3) ! Velocity gradient in VSL
                                                ! (projected along the wall
                                                ! element normal) times the
                                                ! shear viscosity
end type TyParticle_w

! Body elements
type body_element
   integer(4)       :: npart ! Number of body particles
   double precision :: teta_R_IO ! rotation angle for IC and I/O purposes
   double precision :: L_geom(3) ! Element dimensions (Lx,Ly,Lz): a 
                                 ! paralelepiped
   double precision :: dx(3) ! Body particle spacing  
   double precision :: x_CM(3) ! Position of the centre of mass
   double precision :: n_R_IO(3) ! rotation axis for IC and I/O purposes
   integer(4)       :: normal_act(6) ! Boolean value to activate 
                                     ! normals if the element side 
                                     ! is a body side
   double precision :: mass_deact(6) ! (xmin,xmax,ymin,xmax,zmin,zmax) to 
                                     ! deactivate particle masses
end type body_element

! Rigid solid body
type body
   integer(4)       :: npart ! Number of body particles 
   integer(4)       :: Ic_imposed ! Flag to impose Ic in input
   integer(4)       :: n_elem ! Number of body particles 
   integer(4)       :: imposed_kinematics ! Flag for imposed kinematics   
   integer(4)       :: n_records ! Number of records for imposed body 
                                 ! kinematics   
   double precision :: mass ! Mass
   double precision :: umax ! Maximum among the absolute values of the  
                            ! particle velocities
   double precision :: pmax ! Maximum value of pressure
   double precision :: teta_R_IO ! rotation angle for IC and I/O purposes
   double precision :: p_max_limiter ! Maximum pressure value for the maximum 
                                     ! pressure limiter
   double precision :: x_CM(3) ! Position of the centre of mass
   double precision :: alfa(3) ! Rotation angle of the body with respect to 
                               ! the reference system
   double precision :: x_rotC(3) ! Centre of rotation just to configure the
                                 ! initial geometry
   double precision :: u_CM(3) ! Velocity of the centre of mass
   double precision :: omega(3) ! Angular velocity 
   double precision :: Force(3) ! Force
   double precision :: Moment(3) ! Moment / torque
   double precision :: n_R_IO(3) ! rotation axis for IC and I/O purposes
   double precision :: rel_pos_part1_t0(3) ! relative position of the first 
                                           ! particle of the body at the 
                                           ! beginning of the simulation
   double precision :: Ic(3,3) ! Moment of inertia
   double precision :: Ic_inv(3,3) ! Inverse of the moment of inertia
! Array for the imposed body kinematics (n_records*7) 
   double precision,dimension(:,:),allocatable :: body_kinematics 
! Array of the elements of the body 
   type(body_element),dimension(:),allocatable :: elem            
end type body

! Body particle
type body_particle
   integer(4)       :: body ! Body ID 
   integer(4)       :: cell ! Cell ID 
   double precision :: pres ! Pressure
   double precision :: mass ! Mass 
   double precision :: area ! Area (length in 2D) for surface body particles
   double precision :: rel_pos(3) ! Relative position (with respect to the
                                  ! body centre of mass)
   double precision :: pos(3) ! Position (with respect to the origin of the
                              ! reference system)
   double precision :: vel(3) ! Velocity
   double precision :: vel_mir(3) ! Velocity of mirror particles
   double precision :: acc(3) ! Acceleration
   double precision :: normal(3) ! Normal
end type body_particle

! Zone
type TyZone
   logical          :: DBSPH_fictitious_reservoir_flag ! .true.(DB-SPH
                                                       ! fictitious fluid
                                                       ! particles to complete
                                                       ! the kernel support at
                                                       ! free surface in
                                                       ! pre-processing)
   integer(4)       :: ipool 
   integer(4)       :: npoints ! number of topographic vertices used for 
                               ! extrusion (only for a fluid zone extruded from 
                               ! topography)
   integer(4)       :: icol ! Particle colour or number of vertical strips
   integer(4)       :: Medium ! Fluid ID
   integer(4)       :: npointv ! Number of records at imposed kinematics 
                               ! (only for fluid zones)  
   integer(4)       :: IC_source_type ! IC fluid particle distribution from 
                                      ! reservoir vertices and faces (1) or 
                                      ! from Cartesian topography (2)
   integer(4)       :: Car_top_zone ! Zone describing the Cartesian topography,
                                    ! which contains the reservoir 
                                    ! (0 if IC_source_type=1)
   integer(4)       :: plan_reservoir_points ! Number of points describing the
                                             ! reservoir, if(IC_source_type==2)
   integer(4)       :: ID_first_vertex ! ID first vertex of topography 
                                       ! (only for a fluid zone 
                                       ! extruded from topography)
   integer(4)       :: ID_last_vertex ! ID last vertex of topography
                                      ! (only for a fluid zone 
                                      ! extruded from topography)
   integer(4)       :: dam_zone_ID ! ID of the dam zone, related to the
                                   ! eventual reservoir (0 if no dam is 
                                   ! present), if (IC_source_type==2) 
   integer(4)       :: dam_zone_n_vertices ! Number of points describing the
                                           ! horizontal projection of the dam
                                           ! zone, if (dam_zone_ID>1) 
   double precision :: dx_CartTopog ! Spatial resolution of the regular 
                                    ! Cartesian Topography 
   double precision :: H_res ! Height of the reservoir free surface
   double precision :: pool 
   double precision :: trampa ! Time to reach the stationary inlet velocity
   double precision :: valp ! IC for pressure or free surface height
   integer(4)       :: Indix(2) 
   integer(4)       :: limit(2) ! Indices of the first and last particle 
                                ! IDs in the zone 
   double precision :: vel(3) ! Initial velocity 
   double precision :: plan_reservoir_pos(4,2) ! Plane coordinates of the 
                                               ! points (3 or 4), which 
                                               ! describe the reservoir,
                                               ! if (IC_source_type==2)
   double precision :: dam_zone_vertices(4,2) ! Plane coordinates of the points
                                              ! (3 or 4), which describe the 
                                              ! dam zone, if (dam_zone_ID>1)
   double precision :: coordMM(3,2) ! Coordinates of the vertices 
   double precision :: vlaw(0:3,MAXPOINTSVLAW) ! Initial velocity
   character(8)     :: label ! Name 
   character(4)     :: tipo ! Type: "PERI", "SOUR", "OPEN"(, "FLOW", "VELO", 
                            ! "CRIT", "LEVE", "TAPI", "POOL")
   character(1)     :: shape ! Shape 
   character(1)     :: bend ! Colour pattern: uniform or stripped 
   character(2)     :: pressure ! Assigning IC constant pressure or 
                                ! hydrostatic distribution 
   character(3)     :: move ! Motion type: standard or fixed
   character(1)     :: slip ! Boundary slip condition for fixed particles
                            ! (f=free slip; n=no slip)
end type TyZone

! Boundary stretch
type TyBoundaryStretch
   logical          :: laminar_no_slip_check ! Input variable
   integer(4)       :: ColorCode
   integer(4)       :: numvertices
   integer(4)       :: inivertex
   integer(4)       :: iniside
   integer(4)       :: iniface
   integer(4)       :: medium
   integer(4)       :: zone
   double precision :: NormVelocity ! Absolute value of the velocity component,
                                    ! which is normal to the face 
   double precision :: FlowRate ! Flow rate exiting the face 
   double precision :: trampa 
   double precision :: ShearCoeff
   double precision :: velocity(1:SPACEDIM) ! Velocity for "TAPI" zones
   double precision :: PsiCoeff(1:SPACEDIM)
   double precision :: FiCoeff(1:SPACEDIM)
   character(4)     :: tipo ! type: "PERImeter", "SOURce", "OPEN"(, "FIXEd",
                            ! "POOL", "TAPIs", "LEVEl", "FLOW", "VELOcity",
                            ! "CRITic") 
end type TyBoundaryStretch

! Fluid
type TyMedium
   logical :: saturated_medium_flag ! saturated_medium_flag=.true.(fully 
                                    ! saturated medium),.false.(dry medium)  
   integer(4)       :: index ! Fluid ID 
   integer(4)       :: NIterSol ! Number of iterations while the mixture is 
                                ! completely still (beginning of simulation) 
   double precision :: den0 ! IC density
   double precision :: den0_s ! Solid phase reference density 
                              ! (bed-load transport layer)
   double precision :: eps ! Bulk modulus
   double precision :: celerita ! Sound speed 
   double precision :: alfaMon ! Monaghan's alfa parameter 
                               ! (artificial viscosity)
   double precision :: betaMon ! Monaghan's beta parameter 
                               ! (artificial viscosity)
   double precision :: visc ! Kinematic viscosity  
   double precision :: coes ! Cohesion  
   double precision :: limiting_viscosity ! Kinematic viscosity threshold 
                                          ! to save comput. time
   double precision :: numx ! Maximum value for the kinematic viscosity 
                            ! (tuning parameter in Manenti et al., 2012, JHE)
   double precision :: mumx ! Maximum value for the dynamic viscosity 
                            ! (tuning parameter in Manenti et al., 2012, JHE)
   double precision :: taucri ! Critical shear stress (erosion criterion)
   double precision :: cuin ! Exponent of the Constitutive Equation 
   double precision :: phi ! Internal friction angle   
   double precision :: cons ! Consistency   
   double precision :: Cs ! Costant of Smagorinsky's model    
   double precision :: RoughCoef ! Roughness coefficient 
   double precision :: d50 ! 50-th percentile diameter of the granular size
                            ! distribution 
   double precision :: SettlingCoef ! Coefficient for the settling velocity 
                                    ! of the solid grains 
   double precision :: codif ! Diffusion coefficient
   double precision :: Gamma ! Constant in Bachelor equation of state 
                             ! (explosions)
   double precision :: InitialIntEn ! Initial specific internal energy 
   double precision :: gran_vol_frac_max ! Maximum volume fraction of the solid
                                         ! phase within an SPH mixture 
                                         ! particle (bed-load transport layer)
   double precision :: d_90 ! 90-th percentile diameter of the granular size
                            ! distribution 
   character(8)     :: tipo ! Type: "liquid  "(, "gas     ", 
                            ! "general ", "granular", "smagorin")
end type TyMedium

! Boundary side
type TyBoundarySide
   integer(4)       :: stretch
   integer(4)       :: previous_side
   integer(4)       :: vertex(1:SPACEDIM-1)
   integer(4)       :: CloseParticles ! Number of neighbouring particles 
   double precision :: length
   double precision :: CloseParticles_maxQuota ! Maximum height among the  
                                               ! neighbouring particles
   double precision :: angle
   double precision :: velocity(1:SPACEDIM) ! Velocity
   double precision :: T(1:SPACEDIM,1:SPACEDIM)  ! Direction cosines of the
                                                 ! local reference system of
                                                 ! the boundary (the last one
                                                 ! is the normal)   
   double precision :: R(1:SPACEDIM, 1:SPACEDIM)      
   double precision :: RN(1:SPACEDIM, 1:SPACEDIM)
   character(4)     :: tipo ! type: "FIXEd", "PERImeter", "SOURce"(,"TAPIs",  
                            ! "LEVEl", "FLOW", "VELOcity", "CRITic", "OPEN")
end type TyBoundarySide

! Node
type TyNode
   integer(4)       :: name
   double precision :: GX(1:SPACEDIM)
   double precision :: LX(1:SPACEDIM)
end type TyNode

! Face
type TyBoundaryFace
   integer(4)       :: nodes
   integer(4)       :: stretch
   integer(4)       :: CloseParticles ! Number of neighbouring particles
   double precision :: area
   double precision :: CloseParticles_maxQuota ! Maximum height among the  
                                               ! neighbouring particles
   double precision :: T(1:SPACEDIM,1:SPACEDIM) ! Direction cosines of the 
                                                ! local reference system of
                                                ! the boundary (the last 
                                                ! one is the normal)
   double precision :: RPsi(1:SPACEDIM,1:SPACEDIM)
   double precision :: RFi(1:SPACEDIM,1:SPACEDIM)
   double precision :: velocity(1:SPACEDIM) ! Velocity
   type(TyNode)     :: Node(1:MAXFACENODES)
end type TyBoundaryFace

! Close boundary data table
type TyBoundaryData
   integer(4)       :: CloBoNum ! Number of neighbouring faces for a particle
   double precision :: LocXYZ(1:SPACEDIM) ! Face normal vectors in the local 
                                          ! reference system
   double precision :: BoundaryIntegral(1:8) ! Table of integrals  
   double precision :: IntGiWrRdV(1:SPACEDIM,1:SPACEDIM) ! Table of integrals 
end type TyBoundaryData

type TyBoundaryConvexEdge
   double precision :: length                         
   integer(4)       :: face(1:2)
   double precision :: component(1:SPACEDIM)         
   type(TyNode)     :: Node(1:2)
end type TyBoundaryConvexEdge

! Control point
type TyCtlPoint
   integer(4)       :: cella                         
   double precision :: pres ! Interpolated pressure 
   double precision :: dens ! Interpolated density
   double precision :: uni ! Interpolated Shepard coefficient
   double precision :: dist 
   double precision :: coord(3) ! Position
   double precision :: vel(3) ! Interpolated velocity                          
end type TyCtlPoint

! Section
type TySection
   integer(4)       :: ColorCode
   integer(4)       :: Icont(2) ! Pointer to the generated monitoring points
   double precision :: Constant(1:SPACEDIM)          
   double precision :: XYZRange(1:SPACEDIM,2)        
   double precision :: TGLsection(1:SPACEDIM,1:SPACEDIM)  
   character(8)     :: Label                          
   character(2)     :: Tipo ! type: x,y,x,g
end type TySection

! Control line
type TyCtlLine
   integer(4)       :: icont(2) ! Pointers  to the initial and last monitoring  
                                ! points of the discretized line
   character(8)     :: label ! Name 
end type TyCtlLine

! Monitoring section for flow rate (ID depends on the input order)
type tyQ_section_array                                             
   integer(4)       :: n_vertices ! Number of vertices describing a monitoring
                                  ! section for the flow rate (3 or 4)
   double precision :: area ! Area 
   double precision :: normal(3) ! Normal 
   double precision :: vertex(4,3) ! Vertices (in case of 3 vertices, 
                                   ! do not mind about the fourth point)
   double precision :: vertex_loc(4,3) ! Vertices (local coordinates:
                                       ! the origin is vertex1, the third local
                                       ! axis is the normal)
   double precision :: loc_axis(3,3) ! Local axis (axis 1: vertex2-vertex1, 
                                     ! axis 2: vector product of local axis 1
                                     ! and normal; axis 3: normal)
! Section flow rate per fluid type and global: flow_rate(n_fluid_types+1)
   double precision,dimension(:),allocatable  :: flow_rate          
end type

! Monitoring sections for flow rate 
type TyQ_section
   integer(4)       :: n_sect ! Number of monitoring sections for the flow rate 
   integer(4)       :: n_fluid_types ! Number of fluid types (the first ID 
                                     ! fluid types are selected)
   integer(4)       :: it_out_last ! Auxiliary variable to print flow rates 
   double precision :: dt_out ! Writing time step for section flow rates
! type section to monitor flow rates
   type(tyQ_section_array),dimension(:),allocatable :: section           
end type  

! Derived type for bed-load transport layer 
type TyGranular_flows_options
   integer(4)       :: ID_erosion_criterion ! Erosion criterion ID 
                                            ! (1:Shields-Seminara,
                                            ! 2:Shields,3:Mohr-Coulomb)
   integer(4)       :: ID_main_fluid ! ID of the main flow 
   integer(4)       :: monitoring_lines ! Number of monitoring lines aligned 
                                        ! with x- or y-axis
   integer(4)       :: it_out_last ! Auxiliary variable for post-processing
   integer(4)       :: n_max_iterations ! Maximum number of iterations 
                                        ! (iterative process, which estimates
                                        ! u_star)
   integer(4)       :: erosion_flag ! Erosion_flag: 0(activated far from 
                                    ! fronts); 1(not activated), 
                                    ! 2(activated everywhere)
   integer(4)       :: deposition_at_frontiers ! Forced deposition at frontiers
                                               ! (erosion criterion=2): yes(1)
                                               ! or no(2) 
   integer(4)       :: Gamma_slope_flag ! Flag to activate (or not) Gamma_slope       
                                        ! (effects only when ID_erosion 
                                        ! criterion=1) 
   integer(4)       :: saturation_scheme ! saturation_scheme=0(dry soil),1(fully 
                                         ! saturated soil),2(saturation zones 
                                         ! depepending on 
                                         ! time_minimum_saturation and
                                         ! time_maximum_saturation),3(Lagrangian 
                                         ! scheme for saturation conditions)
   double precision :: time_minimum_saturation ! Time related to the minimum 
                                               ! saturation of the granular 
                                               ! material.
   double precision :: time_maximum_saturation ! Time related to the maximum 
                                               ! saturation of the granular 
                                               ! material. So far, 
                                               ! time_minimum_saturation 
                                               ! has to be smaller than (or 
                                               ! equal to) 
                                               ! time_maximum_saturation.
                                               ! When t<=t_min_sat, there is 
                                               ! always phreatic zone below the 
                                               ! free surface and dry soil 
                                               ! elsewhere. When t>=t_max_sat, 
                                               ! the saturation zones are 
                                               ! freezed at t_max_sat.
   double precision :: conv_crit_erosion ! Convergence criterion for erosion  
   double precision :: velocity_fixed_bed ! (optional) velocity_fixed_bed:  
                                          ! velocity threshold (input) to 
                                          ! speed-up a simulation 
   double precision :: dt_out ! Writing time step for the interface monitoring
                              ! lines
   double precision :: x_min_dt ! x_min to involve SPH mixture particles in dt
                                ! estimation
   double precision :: x_max_dt ! x_max to involve SPH mixture particles in dt 
                                ! estimation
   double precision :: y_min_dt ! y_min to involve SPH mixture particles in dt
                                ! estimation
   double precision :: y_max_dt ! y_max to involve SPH mixture particles in dt 
                                ! estimation
   double precision :: z_min_dt ! z_min to involve SPH mixture particles in dt 
                                ! estimation
   double precision :: z_max_dt ! z_max to involve SPH mixture particles in dt 
                                ! estimation 
   double precision :: t_q0 ! t_q0: quake start time 
   double precision :: t_liq ! t_liq: liquefaction time
   logical,dimension(:,:),allocatable :: minimum_saturation_flag ! Free surface 
                                                                 ! flag 
                                                                 ! (presence of 
                                                                 ! the free 
                                                                 ! surface along


                                                                 ! the vertical)
                                                                 ! at the time 
                                                                 ! of minimum 
                                                                 ! saturation of
                                                                 ! the granular 
                                                                 ! material
   logical,dimension(:,:),allocatable :: maximum_saturation_flag ! Free surface 
                                                                 ! flag 
                                                                 ! (presence of 
                                                                 ! the free 
                                                                 ! surface along
                                                                 ! the vertical)
                                                                 ! at the time 
                                                                 ! of maximum 
                                                                 ! saturation of
                                                                 ! the granular 
                                                                 ! material
   integer(4),dimension(:,:),allocatable :: saturation_conditions 
! Saturation conditions: 1 (phreatic zone), 2 (infiltration zone), 3 (dry soil)
! x and/or y coordinates defining the monitoring line/point   
   double precision,dimension(:,:),allocatable :: lines                  
end type

! Face (DB-SPH)
type face_der_type                                                                  
   integer(4),dimension(4)       :: vert_list ! List of face/segment vertices 
                                              ! (triangles in 3D)  
   double precision              :: area ! Face area/length in 3D/2D
   double precision,dimension(3) :: normal ! Face normal 
end type

! Vertex
type vertex_der_type
   double precision,dimension(3) :: pos ! Vertex position                                         
end type

! Surface mesh for DB-SPH boundaries
type DBSPH_surf_mesh_der_type 
! ID of the surface mesh file
   integer(4),allocatable,dimension(:)            :: surface_mesh_file_ID
! List of vertices of the surface mesh
   type(vertex_der_type),allocatable,dimension(:) :: vertices  
! List of faces (both 3D and 2D)
   type(face_der_type),allocatable,dimension(:)   :: faces                    
end type 

! Derived type for DB-SPH bundary treatment scheme
type DBSPH_der_type 
   logical           :: MUSCL_boundary_flag ! Flag to activate (or not) 
                                            ! boundary terms for MUSCL 
                                            ! reconstruction
   logical           :: in_built_monitors ! Flag to activate (or not) 
                                          ! in-built motion of control lines 
   logical           :: Gamma_limiter_flag ! flag to activate or deactivate 
                                           ! Gamma upper limiter (1.)
   logical           :: negative_wall_p_allowed ! pressure of wall elements can 
                                                ! be negative
   logical           :: FS_allowed ! free surface detection can be avoided
   integer(4)        :: n_w ! Number of surface elements
   integer(4)        :: n_monitor_points ! Number of monitoring points
   integer(4)        :: n_monitor_regions ! Number of monitoring regions 
                                          ! (0 or 1) to estimate the Force
                                          ! along x-direction
   integer(4)        :: n_inlet ! Number of inlet sections (to impose DBSPH 
                                ! inlet BC)
   integer(4)        :: n_outlet ! Number of outlet sections (to impose 
                                 ! DBSPH outlet BC)
   integer(4)        :: ply_n_face_vert ! Number of vertices for each surface 
                                        ! mesh face in the .ply input files 
   integer(4)        :: surface_mesh_files ! number of files of the DBSPH 
                                           ! surface meshes
   integer(4)        :: slip_ID ! slip_ID (ID for slip conditions) = 0 
                                ! (free-slip),1 (no-slip), 2 (run-time choice 
                                ! depending on the inner shear viscosity terms 
                                ! in SPH-NS balance equations)
   double precision  :: dx_dxw ! Ratio between the fluid and the 
                               ! semi-particle sizes
   double precision  :: k_w ! Coefficient to compute semi-particle volumes
                            ! (DB-SPH equations)
   double precision  :: monitor_region(6) ! (xmin,xmax,ymin,xmax,zmin,zmax)
                                          ! to detect the monitoring region
! IDs of the monitoring points
   integer(4),allocatable,dimension(:) :: monitor_IDs     
! Number of records for imposed kinematics (surface_mesh_files,n_records)     
   integer(4),allocatable,dimension(:) :: n_kinematics_records 
! rotation_centre(surface_mesh_files,3): centre of rotation for DB-SPH frontiers
   double precision,dimension(:,:),allocatable :: rotation_centre 
! Array for the imposed kinematics of the DBSPH surface elements 
! (surface_mesh_files,n_records,4)
   double precision,dimension(:,:,:),allocatable :: kinematics   
! Array for the inlet sections useful to DB-SPH inlet BC (n_records,10): 
! positions, normal, velocity, length
   double precision,dimension(:,:),allocatable :: inlet_sections 
! Array for the outlet sections useful to DB-SPH outlet BC (n_records,7): 
! positions, normal, length,
   double precision,dimension(:,:),allocatable :: outlet_sections
! DB-SPH surface mesh: vertices and faces  
   type(DBSPH_surf_mesh_der_type) :: surf_mesh                                                                                      
end type 

! Derived type declarations
type(TyGlobal)                 :: Domain
type(TyGriglia)                :: Grid
type(TyParticle)               :: PgZero
type(Tytime_stage)             :: ts_pgZero
type(TyQ_section)              :: Q_sections
type(TyGranular_flows_options) :: Granular_flows_options
type(DBSPH_der_type)           :: DBSPH

end module Hybrid_allocation_module
