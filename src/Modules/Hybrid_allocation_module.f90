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
! Program unit: Hybrid_allocation_module          
! Description: Module to define derived types of both dynamically and statically
!              allocated variables (Di Monaco et al., 2011, EACFM; Manenti et 
!              al., 2012; JHE; Amicarelli et al., 2013, IJNME; Amicarelli et 
!              al., 2015, CAF; Amicarelli et al., 2017, IJCFD).                    
!-------------------------------------------------------------------------------
module Hybrid_allocation_module
use Static_allocation_module
type TyGlobal
   logical :: NormFix
   logical :: Slip
! Output frequency on log file
   integer(4) :: iplot_fr
! Saving frequency                      
   integer(4) :: imemo_fr
! Restart file frequency                      
   integer(4) :: irest_fr
! Monitoring point output frequency                      
   integer(4) :: icpoi_fr
! Output frequency for free surface 
   integer(4) :: ipllb_fr
! Reference fluid for free surface                      
   integer(4) :: ipllb_md
! Initial time step ID
   integer(4) :: istart
! Printing code for log file                        
   integer(4) :: ioutopt
! Maximum step number                       
   integer(4) :: itmax
! Flag for body particles reordering
   integer(4) :: body_part_reorder
#ifdef SPACE_3D
! Rough and slight overestimation of the number of fluid particles in the 
! reservoir (auxiliary parameter useful for extruded reservoir IC)  
   integer(4) :: nag_aux
! Maximum number of neighbouring SA-SPH faces for a computational particle
   integer(4) :: MAXCLOSEBOUNDFACES
! Maximum number of convex edges
   integer(4) :: MAXNUMCONVEXEDGES
#endif
! Density_thresholds flag (default=0; =1 for very low bulk modulus -preliminary 
! simulations-)            
   integer(4) :: density_thresholds
! Flag for Leapfrog scheme
   integer(4) :: time_split
! RK scheme (1,2,3,4) (time_split=0)                    
   integer(4) :: RKscheme
! Stage of RK schemes          
   integer(4) :: time_stage
! Maximum physical time                                           
   double precision :: tmax
! Particle size                           
   double precision :: dx
! h/dx, h=dx*trunc                          
   double precision :: trunc
! dx^D, D: domain dimensionality                       
   double precision :: PVolume                        
   double precision :: coefke  
   double precision :: coefkacl  
   double precision :: CFL
! Viscous stability criterion coefficient, default value: 0.05
   double precision :: vsc_coeff
! Reference pressure
   double precision :: prif
! Output frequency for log file                       
   double precision :: plot_fr
! Frequency for result saving                       
   double precision :: memo_fr
! Output frequency for restart files 
   double precision :: rest_fr
! Output frequency for monitoring points
   double precision :: cpoi_fr
! Output frequency for free surface
   double precision :: pllb_fr 
! Output time step for water depth
   double precision :: depth_dt_out
! Last time step when water depth was printed (auxiliary variable)                         
   double precision :: depth_it_out_last
! Partial smoothing parameter for pressure
   double precision :: TetaP
! Partial smoothing parameter for velocity                         
   double precision :: TetaV
! Kernel support length scale                         
   double precision :: h
! Simulation start time                          
   double precision :: start
! Max number of fluid particles = nag * COEFNMAXPARTI                                  
   double precision :: COEFNMAXPARTI
! maxb(maximum number of fluid neighbours) = COEFNMAXPARTJ * (4h/dx)^D  
   double precision :: COEFNMAXPARTJ
! Time at the beginning of the simulation (origin: beginning of the year)                            
   double precision :: t0
! Time at the beginning of the iterations (origin: beginning of the year)
   double precision :: t_pre_iter
! Gravity acceleration
   double precision :: grav(3)
! Coordinates of 2 vertices of a diagonal of the parallelepiped domain                                            
   double precision :: coord(3,2)
   character(4) :: tipo
   character(len=lencard) :: file
   character(1) :: Psurf
! IC particle distribution noise. "r": slight white noise is added, otherwise 
! no noise.
   character(1) :: RandomPos
end type TyGlobal

! Background positioning grid 
type TyGriglia
! Number of cells
   integer(4) :: nmax
! Number of cells for each axis direction                         
   integer(4) :: ncd(3)
! Grid resolutions                   
   double precision :: dcd(3)
! Coordinates of 2 vertices of a diagonal of the parallelepiped domain
   double precision :: extr(3,2)
end type TyGriglia

! Fluid particle 
type TyParticle
! ID of an inlet neighbouring DB-SPH surface element, if any (otherwise it is
! null)                
   integer(4) :: DBSPH_inlet_ID
! ID of an outlet neighbouring DB-SPH surface element, if any (otherwise it is 
! null)               
   integer(4) :: DBSPH_outlet_ID
! Neighbouring "liquid"/"granular" particle for a computational "granular" / 
! "liquid" particle (erosion criterion)                
   integer(4) :: indneighliqsol
! Index of the mobile ("flu"/mixture) particle, which is the closest to a 
! computational fixed ("sol"/bed) particle and viceversa          
   integer(4) :: ind_neigh_mix_bed
! Bed-load transport flag: =0: generic position; =1: free surface; =2: top of 
! bed-load transport zone; =3: fixed bed; =4: bottom of the granular material; 
! =5: top of the saturated zone.                         
   integer(4) :: blt_flag
! Index of the mobile ("flu"/mixture) particle, which is the closest to the 
! computational mobile granular ("sol"/bed) one     
   integer(4) :: ind_neigh_mob_for_granmob
! Velocity code
   integer(4) :: kodvel
! Density code                    
   integer(4) :: koddens                       
   integer(4) :: CloseBcOut                     
   integer(4) :: cella                          
   integer(4) :: izona
! Colour                         
   integer(4) :: icol
! Fluid ID 
   integer(4) :: imed
! Free Surface (0,3: no free surface; 1,2: free surface; 2,3: deactivated for 
! "Shep" evolution)                           
   integer(4) :: FS
   integer(4) :: laminar_flag
   double precision :: mass
! Density                          
   double precision :: dens
! Pressure                          
   double precision :: pres
! Continuity equation LHS                           
   double precision :: dden
! Shepard coefficient ("Shep": discrete or integral, it depends on the input 
! data)                                                        
   double precision :: uni
! Pressure correction
   double precision :: vpres
! Density correction                          
   double precision :: vden
! Second invariant of the strain rate tensor for incompressible fluids                         
   double precision :: secinv        
   double precision :: dudx
   double precision :: dudy 
   double precision :: dvdx 
   double precision :: dvdy
! (mixture or liquid) kinematic viscosity
   double precision :: kin_visc
! (mixture or liquid) dynamic viscosity
   double precision :: mu
! Stop time                 
   double precision :: tstop                          
   double precision :: mno    
   double precision :: ang
   double precision :: rhoc                           
   double precision :: rhow                          
   double precision :: tiroc                          
   double precision :: cden                           
   double precision :: wden
! Sound speed
   double precision :: Csound
! Density times Shepard coefficient                         
   double precision :: DensShep
! SPH approximation of density at the current time step                     
   double precision :: rhoSPH_new
! SPH approximation of density at the previous time step 
   double precision :: rhoSPH_old
! Lagrangian derivative of Shepard coefficient 
   double precision :: dShep
! Discrete Shepard coefficient
   double precision :: sigma
! Discrete Shepard coefficient involving neighbours of the same fluid   
   double precision :: sigma_same_fluid
! Integral Shepard coefficient
   double precision :: Gamma
! Last value of Gamma before FS=3
   double precision :: Gamma_last_active
! Initial difference between SPH approx. of density and its exact value. This 
! is added in the interior domain to each SPH density approx. to solve inlet 
! problems.
   double precision :: dens_init_err
! Main slope angle of the fixed bed (along the direction aligned with the mean 
! flow; bed-load transport)
   double precision :: Beta_slope
#ifdef SPACE_3D
! Transversal slope angle of the fixed bed (along the direction transversal to 
! the mean flow; bed-load transport)
   double precision :: Gamma_slope
#endif
! Mean of the effective normal stresses
   double precision :: sigma_prime_m
! Pressure of the fluid phase (bed-load transport)
   double precision :: pres_fluid
! Friction velocity representative of the Surface Neutral Boundary Layer 
   double precision :: u_star
#ifdef SPACE_3D
! Lift coefficient for the 3D erosion criterion
   double precision :: C_L
! Drag coefficient for 3D erosion criterion
   double precision :: C_D
#endif
! Ratio (tau/tau_c) between bottom shear stress and critical shear stress 
! (from erosion criterion)
   double precision :: tau_tauc
! Ratio between the 3D critical shear stress (Shields-Seminara) and the 
! analogous 2D criterion value (tauc/tauc,00) 
   double precision :: k_BetaGamma
! Normal of the interface between the mobile particles and the fixed particles 
! (bed-load transport in the presence of an erosion criterion); the normal 
! points inward the mobile domain
   double precision :: normal_int(3)
! Normal of the interface between the mixture and the pure fluid (top of the 
! bed-load transport layer, mixture side)
   double precision :: normal_int_mixture_top(3)
! Normal of the interface between the fully saturated mixture and the rest of 
! the domain
   double precision :: normal_int_sat_top(3)
! Velocity vector before the 3D erosion criterion
   double precision :: vel_old(3)
! normal_int variable before the 3D erosion criterion
   double precision :: normal_int_old(3)
! Density gradient (SPH pseudo-consistent approximation over fluid particles)
   double precision :: drho(3)
! Minimum distance between a SPH mixture particle and a SPH liquid particle
   double precision :: rijtempmin(3)
! Initial velocity for fixed particles
   double precision :: vstart(3)
! Velocity to use for Morris scheme   
   double precision :: velmorr(3)
   double precision :: zer(3)
! Partially smoothed velocity 
   double precision :: var(3)
! Acceleration 
   double precision :: acc(3)
! Position
   double precision :: coord(3)
! Position at the previous time step
   double precision :: CoordOld(3)
! Position at the previous flow rate writing step   
   double precision :: sect_old_pos(3)
! Velocity
   double precision :: vel(3)
! Imposed velocity
   double precision :: velass(3)
! Velocity gradient (SPH pseudo-consistent approximation over fluid particles) 
   double precision :: dvel(3,3)
! BC conditions for fixed particles (f=free-slip; n=no-slip; c=cont slip)
   character(1) :: slip
! Movement type
   character(3) :: vel_type
! Particle "status" ("sol": fixed; "flu": mobile)
   character(3) :: state
end type TyParticle

! Particle intermediate (time stage) values (time integration scheme RK2)
type Tytime_stage
! Stage continuity equation LHS
   double precision :: ts_dden
! Stage density 
   double precision :: ts_dens
! Stage position
   double precision :: ts_coord(3)
! Stage partially smoothed velocity 
   double precision :: ts_var(3)
! Stage velocity
   double precision :: ts_vel(3)
! Stage acceleration
   double precision :: ts_acc(3)
end type Tytime_stage

! Semi-particles and wall elements (DB-SPH)
type TyParticle_w
! Cell number
   integer(4) :: cella
! Index of the adjacent_faces (3 at maximum)
   integer(4) :: adjacent_faces(3)
! wet=1 for wet wall element (distance from fluid particle smaller than 1.3dx)
   integer(4) :: wet
! ID of the surface mesh file
   integer(4) :: surface_mesh_file_ID
! Density
   double precision :: dens
! Pressure
   double precision :: pres
! Area(3D)/length(2D) of the wall element
   double precision :: weight
! Mass of the semi-particle 
   double precision :: mass
! Depth coefficient
   double precision :: k_d
! Semi-particle volume (area in 2D)
   double precision :: volume
! Discrete Shepard coefficient of the wall elements depending on fluid 
! particles (not on semi-particles) 
   double precision :: sigma
! Kinematic viscosity of the semi-particle
   double precision :: kin_visc_semi_part                             
! Normal
   double precision :: normal(3)
! Position
   double precision :: coord(3)
! Velocity
   double precision :: vel(3)
! Velocity gradient in VSL (projected along the wall element normal) times the
! shear viscosity
   double precision :: grad_vel_VSL_times_mu(3)
end type TyParticle_w

! Body elements
type body_element
! Number of body particles
   integer(4) :: npart
! rotation angle for IC and I/O purposes
   double precision :: teta_R_IO
! Element (parallelepiped) dimensions (Lx,Ly,Lz)
   double precision :: L_geom(3)
! Body particle spacing
   double precision :: dx(3)
! Position of the centre of mass
   double precision :: x_CM(3)
! rotation axis for IC and I/O purposes
   double precision :: n_R_IO(3)
! Boolean value to activate normals if the element side is a body side
   integer(4) :: normal_act(6)
! (xmin,xmax,ymin,xmax,zmin,zmax) to deactivate particle masses
   double precision :: mass_deact(6)
end type body_element

! Rigid solid body
type body
! Number of body particles
   integer(4) :: npart
! Flag to impose Ic in input
   integer(4) :: Ic_imposed
! Number of body particles
   integer(4) :: n_elem
! Flag for imposed kinematics
   integer(4) :: imposed_kinematics
! Number of records for imposed body kinematics      
   integer(4) :: n_records
! Mass
   double precision :: mass
! Maximum among the absolute values of the particle velocities
   double precision :: umax
! Maximum value of pressure
   double precision :: pmax
! rotation angle for IC and I/O purposes
   double precision :: teta_R_IO
! Maximum pressure value for the maximum pressure limiter
   double precision :: p_max_limiter
! Position of the centre of mass
   double precision :: x_CM(3)
! Rotation angle of the body with respect to the reference system
   double precision :: alfa(3)
! Centre of rotation just to configure the initial geometry
   double precision :: x_rotC(3)
! Velocity of the centre of mass
   double precision :: u_CM(3)
! Angular velocity
   double precision :: omega(3)
! Force 
   double precision :: Force(3)
! Moment / torque
   double precision :: Moment(3)
! Rotation axis for IC and I/O purposes
   double precision :: n_R_IO(3)
! Relative position of the first particle of the body at the beginning of the 
! simulation
   double precision :: rel_pos_part1_t0(3)
! Moment of inertia
   double precision :: Ic(3,3)
! Inverse of the moment of inertia
   double precision :: Ic_inv(3,3)
! Array for the imposed body kinematics (n_records*7) 
   double precision,dimension(:,:),allocatable :: body_kinematics 
! Array of the elements of the body 
   type(body_element),dimension(:),allocatable :: elem            
end type body

! Body particle
type body_particle
! Body ID
   integer(4) :: body
! Cell ID
   integer(4) :: cell
! Pressure
   double precision :: pres
! Mass
   double precision :: mass
! Area (length in 2D) for surface body particles
   double precision :: area
! Relative position (with respect to the body barycentre)
   double precision :: rel_pos(3)
! Position (with respect to the origin of the reference system)
   double precision :: pos(3)
! Velocity
   double precision :: vel(3)
! Velocity of mirror particles
   double precision :: vel_mir(3)
! Acceleration
   double precision :: acc(3)
! Normal
   double precision :: normal(3)
end type body_particle

! Zone
type TyZone
! Flag to activate/deactivate DB-SPH fictitious fluid particles to complete the 
! kernel support at free surface to impose initial conditions
   logical :: DBSPH_fictitious_reservoir_flag
   integer(4) :: ipool
! Number of topographic vertices used for extrusion (only for a fluid zone 
! extruded from topography)
   integer(4) :: npoints
! Particle colour or number of vertical strips
   integer(4) :: icol
! Fluid ID
   integer(4) :: Medium
! Number of records for kinematics imposed (only for fluid zones)  
   integer(4) :: npointv
! IC fluid particle distribution from reservoir vertices and faces 
! (IC_source_type=1) or from Cartesian topography (IC_source_type=2)
   integer(4) :: IC_source_type
! Zone describing the Cartesian topography, which contains the reservoir 
! (Car_top_zone=0 when IC_source_type=1)
   integer(4) :: Car_top_zone
   integer(4) :: slip_coefficient_mode
#ifdef SPACE_3D   
! Number of points describing the reservoir if (IC_source_type==2)
   integer(4) :: plan_reservoir_points
! ID of the first vertex of topography (for a fluid zone extruded from 
! topography)
   integer(4) :: ID_first_vertex
! ID of the last vertex of topography (for a fluid zone extruded from 
! topography) 
   integer(4) :: ID_last_vertex
! ID of the dam zone, related to the eventual reservoir (dam_zone_ID=0 if no 
! dam is present), if(IC_source_type==2) 
   integer(4) :: dam_zone_ID
! Number of points describing the horizontal projection of the dam zone, 
! if (dam_zone_ID>1) 
   integer(4) :: dam_zone_n_vertices
! Spatial resolution of the regular Cartesian Topography 
   double precision :: dx_CartTopog
! Height of the reservoir free surface
   double precision :: H_res
#endif
! Input variable for the boundary shear stress: 
!    slip coefficient (if slip_coefficient_mode==1)
!    mean diameter of the wall roughness elements (if slip_coefficient_mode==2) 
   double precision :: BC_shear_stress_input
! Average computed slip coefficient
   double precision :: avg_comp_slip_coeff
   double precision :: pool
! IC for pressure or free surface height 
   double precision :: valp
   integer(4) :: Indix(2)
! Indices of the first and last particle IDs in the zone  
   integer(4) :: limit(2)
! Initial velocity
   double precision :: vel(3)
! Coordinates of the vertices
   double precision :: coordMM(3,2)
! Horizontal coordinates of the points (3 or 4), which describe the reservoir,
! if (IC_source_type==2)
   double precision :: plan_reservoir_pos(4,2)
#ifdef SPACE_3D
! Horizontal coordinates of the points (3 or 4), which describe the dam zone, 
! if (dam_zone_ID>1)
   double precision :: dam_zone_vertices(4,2)
#endif
! Initial velocity
   double precision :: vlaw(0:3,MAXPOINTSVLAW)
   character(1) :: shape
! Colour pattern: uniform or stripped 
   character(1) :: bend
! Assigning IC constant pressure or hydrostatic distribution 
   character(2) :: pressure
! Motion type: standard or fixed
   character(3) :: move
! Type: "PERI", "SOUR", "OPEN"(, "FLOW", "VELO", "CRIT", "LEVE", "TAPI", "POOL") 
   character(4) :: tipo
! Name
   character(8) :: label
end type TyZone

! Boundary stretch
type TyBoundaryStretch
   logical :: laminar_no_slip_check
   integer(4) :: ColorCode
   integer(4) :: numvertices
   integer(4) :: inivertex
#ifdef SPACE_3D
   integer(4) :: iniface
#elif defined SPACE_2D
   integer(4) :: iniside
#endif
   integer(4) :: medium
   integer(4) :: zone
! Absolute value of the velocity component which is normal to the face 
   double precision :: NormVelocity
! Flow rate exiting the face
   double precision :: FlowRate
! Velocity for "TAPI" zones
   double precision :: velocity(1:SPACEDIM)
   double precision :: PsiCoeff(1:SPACEDIM)
   double precision :: FiCoeff(1:SPACEDIM)
! Type: "PERI", "SOUR", "OPEN"(, "FLOW", "VELO", "CRIT", "LEVE", "TAPI", "POOL")
   character(4) :: tipo
end type TyBoundaryStretch

! Fluid
type TyMedium
! saturated_medium_flag=.true.(fully saturated medium),.false.(dry medium)
   logical :: saturated_medium_flag
! Fluid ID
   integer(4) :: index
! Number of iterations while the mixture is completely still (IC) 
   integer(4) :: NIterSol
! Fluid reference density
   double precision :: den0
! Solid phase reference density
   double precision :: den0_s
! Bulk modulus
   double precision :: eps
! Sound speed
   double precision :: celerita
! Monaghan's alfa parameter (artificial viscosity) 
   double precision :: alfaMon
! Monaghan's beta parameter (artificial viscosity)
   double precision :: betaMon
! Kinematic viscosity (in case of granular media, kin_visc=mumx)
   double precision :: kin_visc
   double precision :: limiting_viscosity
! Maximum value for the dynamic viscosity to detect the elasto-plastic regime
   double precision :: mumx
! Internal friction angle
   double precision :: phi
! Roughness coefficient for the erosion criteria
   double precision :: RoughCoef
! 50-th percentile diameter of the granular size distribution 
   double precision :: d50
! Maximum volume fraction of the solid phase within an SPH mixture particle
   double precision :: gran_vol_frac_max
! 90-th percentile diameter of the granular size distribution
   double precision :: d_90
! Type: "liquid  ","granular" 
   character(8) :: tipo
end type TyMedium

#ifdef SPACE_2D
! Boundary side
type TyBoundarySide
   integer(4) :: stretch
   integer(4) :: previous_side
   integer(4) :: vertex(1:SPACEDIM-1)
! Number of neighbouring particles
   integer(4) :: CloseParticles 
   double precision :: length
! Maximum height among the neighbouring particles
   double precision :: CloseParticles_maxQuota
   double precision :: angle
   double precision :: velocity(1:SPACEDIM)
! Direction cosines of the local reference system of the boundary 
! (T(SPACEDIM,1:SPACEDIM) is the normal) 
   double precision :: T(1:SPACEDIM,1:SPACEDIM)
   double precision :: R(1:SPACEDIM,1:SPACEDIM)      
   double precision :: RN(1:SPACEDIM,1:SPACEDIM)
! Type: "PERI", "SOUR", "OPEN"(, "FLOW", "VELO", "CRIT", "LEVE", "TAPI", "POOL")
   character(4) :: tipo
end type TyBoundarySide
#endif

! Node
type TyNode
   integer(4) :: name
! Position of the vertex in the global reference system 
   double precision :: GX(1:SPACEDIM)
! Position of the vertex in the local reference system
   double precision :: LX(1:SPACEDIM)
end type TyNode

#ifdef SPACE_3D
! Face
type TyBoundaryFace
   integer(4) :: nodes
   integer(4) :: stretch
! Number of neighbouring particles
   integer(4) :: CloseParticles
   double precision :: area
! Maximum height among the neighbouring particles
   double precision :: CloseParticles_maxQuota
! Direction cosines of the local reference system of the boundary 
! (T(SPACEDIM,1:SPACEDIM) is the normal)
   double precision :: T(1:SPACEDIM,1:SPACEDIM)
   double precision :: RPsi(1:SPACEDIM,1:SPACEDIM)
   double precision :: RFi(1:SPACEDIM,1:SPACEDIM)
   double precision :: velocity(1:SPACEDIM)
   type(TyNode) :: Node(1:MAXFACENODES)
end type TyBoundaryFace
#endif

! Close boundary data table
type TyBoundaryData
! Number of neighbouring faces for a particle
   integer(4) :: CloBoNum
! Neighbouring particle coordinates in the frontier local reference system
   double precision :: LocXYZ(1:SPACEDIM)
! Table of integrals
   double precision :: BoundaryIntegral(1:8)
! Table of integrals
   double precision :: IntGiWrRdV(1:SPACEDIM,1:SPACEDIM)
end type TyBoundaryData

#ifdef SPACE_3D
type TyBoundaryConvexEdge
   double precision :: length                         
   integer(4) :: face(1:2)
   double precision :: component(1:SPACEDIM)         
   type(TyNode) :: Node(1:2)
end type TyBoundaryConvexEdge
#endif

! Control point
type TyCtlPoint
   integer(4) :: cella
! Interpolated pressure                    
   double precision :: pres
! Interpolated density 
   double precision :: dens
! Interpolated Shepard coefficient
   double precision :: uni
   double precision :: dist
! Position
   double precision :: coord(3)
! Interpolated velocity
   double precision :: vel(3)
end type TyCtlPoint

! Control line
type TyCtlLine
! Pointers  to the initial and last monitoring points of the discretized line
   integer(4) :: icont(2)
   character(8) :: label
end type TyCtlLine

#ifdef SPACE_3D
! Monitoring section for flow rate (ID depends on the input order)
type tyQ_section_array
! Number of vertices describing a monitoring section for the flow rate (3 or 4)
   integer(4) :: n_vertices
   double precision :: area
   double precision :: normal(3)
! Vertices (in case of 3 vertices, do not mind about the fourth point)
   double precision :: vertex(4,3)
! Vertices (local coordinates: the origin is vertex1, the third local axis is 
! the normal)
   double precision :: vertex_loc(4,3)
! Local axis (axis 1: vertex2-vertex1, axis 2: vector product of local axis 1
! and the normal; axis 3: normal)
   double precision :: loc_axis(3,3)
! Section flow rate per fluid type and global: flow_rate(n_fluid_types+1)
   double precision,dimension(:),allocatable  :: flow_rate          
end type

! Monitoring sections for flow rate 
type TyQ_section
! Number of monitoring sections for the flow rate
   integer(4) :: n_sect
! Number of fluid types (the first ID fluid types are selected)
   integer(4) :: n_fluid_types
! Auxiliary variable to print flow rates
   integer(4) :: it_out_last
! Writing time step for the monitoring sections
   double precision :: dt_out
   type(tyQ_section_array),dimension(:),allocatable :: section           
end type

! Substation array
type type_substation
   integer(4) :: n_vertices
! Substation type
   integer(4) :: type_ID
! Number of DEM vertices within the substation surface
   integer(4) :: n_DEM_vertices
   double precision :: area
! Sum (over time) of the filtered values of POS (useful to estimate EOS); two 
! values are available depending on the assessment method: one assumes that 
! Ysub alternatively depends on two estimations of the water depth. The first 
! is a raw estimation; the latter filters the atomization and the wave breaking 
! effects.
   double precision :: POS_fsum(2)
! Maximum substation fluid/mixture depth; two values are available depending on 
! the assessment method. 
   double precision :: Ymax(2)
! Expected Outage Time; two values are available depending on the assessment 
! method.
   double precision :: EOT(2)
! Substation value (euros)
   double precision :: Val
! List of DEM vertices associated with the substation
   integer(4),dimension(:),allocatable :: DEMvert
! Vertices of the polygon which describes the substation
   double precision :: vert(6,2)
end type

! Monitoring the electrical substations (ref. input file template)
type type_substations
   integer(4) :: n_sub
! Last time step of substation result writing (auxiliary variable)
   integer(4) :: it_out_last
! Writing time step for substations
   double precision :: dt_out
! type for the substation array
   type(type_substation),dimension(:),allocatable :: sub          
end type
#endif

! Derived type for bed-load transport layer 
type TyGranular_flows_options
! Dense granular flow configuration
   integer(4) :: KTGF_config
   integer(4) :: ID_main_fluid
! Number of monitoring lines aligned with x- or y-axis
   integer(4) :: monitoring_lines
! Auxiliary variable for post-processing
   integer(4) :: it_out_last
! Maximum number of iterations (iterative process, which estimates u_star)
   integer(4) :: n_max_iterations
! Erosion_flag: 0(activated far from fronts); 1(not activated), 2(activated 
! everywhere)
   integer(4) :: erosion_flag
! Forced deposition at frontiers 
   integer(4) :: deposition_at_frontiers
#ifdef SPACE_3D
! Flag to activate (or not) Gamma_slope
   integer(4) :: Gamma_slope_flag
#endif
   integer(4) :: saturation_scheme
   double precision :: time_minimum_saturation
   double precision :: time_maximum_saturation
! Convergence criterion for erosion
   double precision :: conv_crit_erosion
! Velocity threshold (input)
   double precision :: velocity_fixed_bed
   double precision :: dt_out
! x_min to involve SPH mixture particles in dt estimation
   double precision :: x_min_dt
   double precision :: x_max_dt
   double precision :: y_min_dt
   double precision :: y_max_dt
   double precision :: z_min_dt
   double precision :: z_max_dt
! Quake start time
   double precision :: t_q0
! Liquefaction time
   double precision :: t_liq
! Free surface flag at the time of minimum saturation of the granular material
   logical,dimension(:,:),allocatable :: minimum_saturation_flag
   logical,dimension(:,:),allocatable :: maximum_saturation_flag
   integer(4),dimension(:,:),allocatable :: saturation_conditions
! x and/or y coordinates defining the monitoring line/point   
   double precision,dimension(:,:),allocatable :: lines                  
end type

! Face (DB-SPH)
type face_der_type
! List of face/segment vertices (triangles in 3D)                                             
   integer(4),dimension(4) :: vert_list  
 ! Face area/length in 3D/2D
   double precision :: area
   double precision,dimension(3) :: normal 
end type

! Vertex
type vertex_der_type
! Position
   double precision,dimension(3) :: pos                                         
end type

! Surface mesh for DB-SPH boundaries
type DBSPH_surf_mesh_der_type 
! ID of the surface mesh file
   integer(4),allocatable,dimension(:) :: surface_mesh_file_ID
! List of vertices of the surface mesh
   type(vertex_der_type),allocatable,dimension(:) :: vertices  
! List of faces (both 3D and 2D)
   type(face_der_type),allocatable,dimension(:) :: faces                    
end type 

! Derived type for DB-SPH bundary treatment scheme
type DBSPH_der_type
! Flag to activate (or not) boundary terms for MUSCL reconstruction
   logical :: MUSCL_boundary_flag
! Flag to activate (or not) in-built motion of control lines 
   logical :: in_built_monitors
! Flag to activate or deactivate Gamma upper limiter (1.)
   logical :: Gamma_limiter_flag
! Pressure of wall elements can be negative
   logical :: negative_wall_p_allowed
! Free surface detection can be avoided
   logical :: FS_allowed
! Number of surface elements
   integer(4) :: n_w
! Number of monitoring points
   integer(4) :: n_monitor_points
! Number of monitoring regions (0 or 1) to estimate the Force along x-direction
   integer(4) :: n_monitor_regions
! Number of inlet sections (to impose DBSPH inlet BC)
   integer(4) :: n_inlet
! Number of outlet sections (to impose DBSPH outlet BC)
   integer(4) :: n_outlet
! Number of vertices for each surface mesh face in the ".ply" input files 
   integer(4) :: ply_n_face_vert
! Number of files of the DBSPH surface meshes
   integer(4) :: surface_mesh_files
! slip_ID (ID for slip conditions) = 0 (free-slip),1 (no-slip), 2 (run-time 
! choice depending on the inner shear viscosity terms in SPH-NS balance 
! equations)  
   integer(4) :: slip_ID
! Ratio between the fluid and the semi-particle sizes
   double precision :: dx_dxw
! Coefficient to compute semi-particle volumes (DB-SPH equations)
   double precision :: k_w
! (xmin,xmax,ymin,xmax,zmin,zmax) to detect the monitoring region
   double precision :: monitor_region(6)
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
type(TyGlobal) :: Domain
type(TyGriglia) :: Grid
type(TyParticle) :: PgZero
type(Tytime_stage) :: ts_pgZero
#ifdef SPACE_3D
type(TyQ_section) :: Q_sections
type(type_substations) :: substations
#endif
type(TyGranular_flows_options) :: Granular_flows_options
type(DBSPH_der_type) :: DBSPH

end module Hybrid_allocation_module
