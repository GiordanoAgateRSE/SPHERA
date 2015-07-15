
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : AdM_User_Type
!
! Last updating : April 18, 2013
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Add the message unit on the screen
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  30Jun11        Time stage parameters
! 04  Amicarelli/Agate  30Nov11        BSPH: wall element parameters
! 05  Amicarelli/Agate  13nov12        (AA501b) Body dynamics
! 06  Amicarelli/Agate  18apr12        add flag for wall elements reordering and maximum number of neighbours read in input
! 07  Amicarelli        8apr14         (v5.04) New variables for: water depth, memory management from input data, real-time forecast of the simulation duration, 
!                                      3D erosion criterion (also with mixture - fixed bed interactions), granular flows, IC reservoir + dam in input from PV, hydrographs 
! 08  Amicarelli        26Jan15        DBSPH-input (AA601) . Derived types for DBSPH: vertex, face, surf_mesh, DBSPH. Modificatuions of the drived type Typarticle_w.
!
!************************************************************************************
! Module purpose : Global variables and data structure declaration
!
! Calling routine: all
!
! Called routines: none
!
!************************************************************************************
!
module AdM_User_Type
!
!------------------------------------------------------------------------------
! additional user-defined types specific for semianalitical boundary conditions
!------------------------------------------------------------------------------
!
use GLOBAL_MODULE
!
!
!.. defines the structure of global variables
!
type TyGlobal
   character(4)           :: tipo
   character(80)          :: file
   character(1)           :: Psurf                         ! Smoothing pressione sulla superficie  O=original, M=deltaP con Pidrostatica, A=pesato con Patmosfera
   character(1)           :: RandomPos                     ! Posizionamento iniziale delle particelle  niente=originale R=Random
   integer(4)             :: iplot_fr                      ! frequenza di scrittura dei risultati su file out
   integer(4)             :: imemo_fr                      ! frequenza di memorizzazione dei risultati
   integer(4)             :: irest_fr                      ! frequenza di scrittura del file di restart
   integer(4)             :: icpoi_fr                      ! frequenza di memorizzazione quantita' nei punti sonda
   integer(4)             :: ipllb_fr                      ! frequenza di memorizzazione pelo libero
   integer(4)             :: ipllb_md                      ! mezzo di riferimento per il calcolo del pelo libero
   integer(4)             :: istart                        ! passo iniziale della simulazione
   integer(4)             :: ioutopt                       ! codice per stampa su file di output
   integer(4)             :: itmax                         ! numero di step massimo della simulazione
   double precision       :: coord(3,2)                    ! coordinate degli estremi della diagonale
   double precision       :: tmax                          ! tempo massimo della simulazione
   double precision       :: dd                            ! valore della scala di discretizzazione
   double precision       :: trunc                         ! 2h=dd*trunc
   double precision       :: PVolume                       ! Volume della particella =dd^2 in 2D o dd^3 in 3D
   double precision       :: coefke                        ! 
   double precision       :: coefkacl                      ! 
!AA401 start
!   double precision       :: cote                          ! parametro di riduzione del passo di tempo
   double precision       :: CFL                           ! CFL number
   integer                :: time_split                    ! option to activate the staggered Euler time integration
   integer                :: RKscheme                      ! RK scheme (1,2,3,4); valid if time_split=0
!AA401 end
!AA402 start
   integer                :: time_stage                    ! stage of RK schemes
!AA402 end
   double precision       :: grav(3)                       ! componenti dell'accelerazione di gravita'
   double precision       :: prif                          ! pressione di riferimento del sistema
   double precision       :: plot_fr                       ! intervallo di tempo per la scrittura dei risultati su file out
   double precision       :: memo_fr                       ! intervallo di tempo per la memorizzazione dei risultati
   double precision       :: rest_fr                       ! intervallo di tempo per scrittura file di restart
   double precision       :: cpoi_fr                       ! intervallo di tempo per memo quantita' nei punti sonda
   double precision       :: pllb_fr                       ! intervallo di tempo per memo quantita' per il pelo libero
!AA504_B4 start
   double precision       :: depth_dt_out                  ! Writing time step (s) for water depth
   double precision       :: depth_it_out_last             ! Last time step when water depth was printed (auxiliary variable)
!AA504_B4 end   
   double precision       :: TetaP                         ! valore del parametro di smoothing sulla pressione
   double precision       :: TetaV                         ! valore del parametro di smoothing sulla velocita'
   double precision       :: pre                           ! valore della pressione di riferimento
   double precision       :: h                             ! valore del semiraggio di influenza delle particelle
   double precision       :: start                         ! istante di inizio della simulazione
   logical                :: NormFix                       ! 
   logical                :: Slip                          ! 
!AA504
   double precision       :: COEFNMAXPARTI                 ! COEFNMAXPARTI: max_neumber of fluid particles=nag*COEFNMAXPARTI
!AA503 start
!AA504sub comm
   double precision       :: COEFNMAXPARTJ                 ! COEFNMAXPARTJ: maxb(maximum number of fluid neighbours) =  COEFNMAXPARTJ * (4h/dx)^D(domain dimension)
   integer(4)             :: body_part_reorder             ! flag for body particles reordering 
   integer(4)             :: MAXCLOSEBOUNDFACES            ! Maximum number of neighbouring SASPH faces for a computational particle
!AA504
   integer(4)             :: MAXNUMCONVEXEDGES             ! Maximum number of convex edges 
   integer(4)             :: density_thresholds            ! density_thresholds flag (default=0; =1 when E is chosen very low for preliminary simulations)
!AA503 end
!AA504 start
   double precision       :: t0                            ! time (s) at the beginning of the simulation (origin: beginning of the year) 
   double precision       :: t_pre_iter                    ! time (s) at the beginning of the iterations (origin: beginning of the year) 
!AA504 end
end type TyGlobal
!
!
!.. defines the structure of the virtual computational grid
!
type TyGriglia
   integer(4)             :: ncd(3)                        ! numero di celle nelle tre direzioni
   integer(4)             :: nmax                          ! numero totale di celle della griglia
   double precision       :: dcd(3)                        ! deltaX, deltaY, deltaZ delle celle
   double precision       :: extr(3,2)                     ! coordinate dei due estremi della griglia
end type TyGriglia
!
!
!.. defines the structure of the particle values and attributes
!
type TyParticle
   character(3)           :: vel_type                      ! tipo di moto
   character(1)           :: slip                          ! condizione sullo scorrimento (solo particelle fisse) f=free slip;  n=no slip;  c=cont slip
   character(3)           :: state                         ! stato della particella 'sol' solido (non si muove) 'flu' fluido (si muove)
!AA601 start
   integer(4)             :: DBSPH_inlet_ID                ! ID of a inlet neighbouring DBSPH surface element if any, otherwise it is null 
   integer(4)             :: DBSPH_outlet_ID               ! ID of a outlet neighbouring DBSPH surface element if any, otherwise it is null
!AA601 end
   integer(4)             :: indneighliqsol                ! indice particella 'liquid' vicino ad una 'granular' o viceversa
!AA504 start
   integer(4)             :: ind_neigh_mix_bed             ! index of the mobile ("flu"/mixture) particle which is the closest to a computational fixed ("sol"/bed) particle 
                                                           ! and viceversa
   integer(4)             :: blt_flag                      ! bed load transport flag: 0:generic position, 1:free surface; 2:top of bed load transport zone, 3:fixed bed 
   integer(4)             :: ind_neigh_mob_for_granmob     ! index of the mobile ("flu"/mixture) particle which is the closest to the computational mobile granular ("sol"/bed) one
!AA504 end
   integer(4)             :: kodvel                        ! codice velocita'
   integer(4)             :: koddens                       ! codice densita'
!!!   integer(4)             :: npar2h                        ! numero di particelle entro 2h  ??? non utilizzata ???
   integer(4)             :: CloseBcOut                    ! condizione al contorno vicina di uscita
   integer(4)             :: cella                         ! numero della cella a cui appartiene la particella
   integer(4)             :: izona                         ! zona di appartenenza della particella
   integer(4)             :: icol                          ! colore della particella
   integer(4)             :: imed                          ! indice del mezzo di appartenenza
   double precision       :: coord(3)                      ! posizione 
   double precision       :: CoordOld(3)                   ! posizione al passo precedente
!AA504
   double precision       :: sect_old_pos(3)               ! coordinates of fluid particles at the previous flow rate writing step   
   double precision       :: vel(3)                        ! velocita'
   double precision       :: velass(3)                     ! velocita' assegnata
   double precision       :: densass                       ! densita' assegnata
   double precision       :: acc(3)                        ! accelerazione
   double precision       :: mass                          ! massa
   double precision       :: dens                          ! densita'
   double precision       :: pres                          ! pressione
   double precision       :: dden                          ! variazione di densita'
   double precision       :: uni                           ! stima dell'unita' !AA406: Shepard's coefficient (Shep: discrete or integral, it depends on the case)
   double precision       :: zer(3)                        ! normale alla particella (solo per particelle fisse)
   double precision       :: var(3)                        ! velocita' corretta
   double precision       :: vpres                         ! correzione di pressione
   double precision       :: vden                          ! correzione di densita'
!AA504 sub comm   
   double precision       :: secinv                        ! Second invariant of the strain rate tensor for incompressible fluids 
   double precision       :: dudx                          ! derivata della comp. oriz. della velocita' rispetto ad x
   double precision       :: dudy                          ! derivata della comp. oriz. della velocita' rispetto ad y
   double precision       :: dvdx                          ! derivata della comp. vert. della velocita' rispetto ad x
   double precision       :: dvdy                          ! derivata della comp. vert. della velocita' rispetto ad y
   double precision       :: visc                          ! viscosita' cinematica
!AA504 sub comm   
   double precision       :: mu                            ! dynamic viscosity; equivalent mixture dynamic viscosity in case of granular flows
   double precision       :: vstart(3)                     ! velocita' iniziale per particelle fisse
   double precision       :: velmorr(3)                    ! velocita' da utilizzare per lo schema di Morris
   double precision       :: tstop                         ! tempo di arresto
   double precision       :: mno                           ! 
   double precision       :: ang                           ! angolo
!AA504 sub   
   double precision       :: rijtempmin(3)                 ! minimum distance between a SPH mixture particle and a SPH pure fluid particle
   !!!! parti modello bifluido
   double precision       :: VolFra                        ! Volume Fraction
   double precision       :: rhoc                          ! 
   double precision       :: rhow                          ! 
   double precision       :: tiroc                         ! rho tilde nel granulare diffusione
   double precision       :: cden                          ! 
   double precision       :: wden                          ! 
   double precision       :: diffu                         ! 
   double precision       :: coefdif                       ! Coefficiente diffusivo
   double precision       :: veldif(3)                     ! velocita' modello diffusione
   !!!! parti modello bifluido
   !! modello esplosioni
   double precision       :: IntEn                         ! Specific Internal Energy 
   double precision       :: Envar                         ! Specific Internal Energy sums contributions for smoothing
   double precision       :: dEdT                          ! variazione Specific Internal Energy
   double precision       :: Csound                        ! Sound speed
!
!AA406 start
   double precision       :: drho(3)                       ! density gradient (SPH pseudo-consistent approx. over fluid part.)
   double precision       :: dvel(3,3)                     ! velocity gradient (SPH pseudo-consistent approx. over fluid part.)
!AA406 end
!
!AA406test
   double precision       :: DensShep                      ! Density x Shepard's coefficient
   double precision       :: rhoSPH_new                    ! SPH approximation of density at the new time step 
   double precision       :: rhoSPH_old                    ! SPH approximation of density at the old time step 
   double precision       :: dShep                         ! Lagrangian derivative of the Shepard coefficient Shep or uni
   double precision       :: sigma                         ! Discrete Shepard coefficient
   double precision       :: Gamma                         ! Integral Shepard coefficient
   integer(4)             :: FS                            ! Free Surface (0,3: no free surface; 1,2: free surface)
                                                           !              (2,3: deactivated for Shep evolution)
   double precision       :: Gamma_last_active             ! last value of Gamma before FS=3
   double precision       :: dens_init_err                 ! initial difference between SPH approx of density and its exact value
                                                           ! this is added in the interior domain to each SPH density approx 
                                                           ! to solve inlet problems
!AA504 start
   double precision       :: Beta_slope                    ! main slope angle of the fixed bed (along the direction aligned with the mean flow; granular flows)
   double precision       :: Gamma_slope                   ! transversal slope angle of the fixed bed (along the direction transversal to the mean flow; granular flows)
   double precision       :: sigma_prime                   ! effective normal stress (granular flows)
   double precision       :: u_star                        ! Friction velocity representative of the Surface Neutral Boundary Layer (granular flows) 
   double precision       :: C_L                           ! Lift coefficient for 3D erosion criterion (granular flows)
   double precision       :: C_D                           ! Drag coefficient for 3D erosion criterion (granular flows)
   double precision       :: tau_tauc                      ! Ratio (tau/tau_c) between bottom shear stress and critical shear stress (from erosion criterion)
   double precision       :: k_BetaGamma                   ! Ratio between 3D critical shear stress (Shields-Seminara) and analogous 2D criterion (tauc/tauc,00) 
   double precision       :: Bn                            ! Bingham number   
   double precision       :: normal_int(3)                 ! Normal of the interface between the granular mixture and the pure fluid (pointing outward the main flow domain)
   double precision       :: vel_old(3)                    ! Velocity vector before the erosion criterion (3D erosion criterion)
   double precision       :: normal_int_old(3)             ! (fixed bed - mixture) interface normal before the erosion criterion (3D erosion criterion)
!A504 end
!! modello esplosioni
!25gen2011   logical                :: punta                         ! riconosce le sporgenze solide
end type TyParticle
!
!AA402
!.. defines the structure of the particle intermediate (time stage) values (time integrations RK2/3/4)
!
type Tytime_stage
   double precision       :: ts_vel(3)                     ! stage velocity
   double precision       :: ts_acc(3)                     ! stage acceleration
   double precision       :: ts_dden                       ! stage continuity equation right hand side 
   double precision       :: ts_dens                       ! stage density
   double precision       :: ts_var(3)                     ! stage corrected velocity
   double precision       :: ts_IntEn                      ! stage specific internal energy 
   double precision       :: ts_dEdT                       ! stage energy equation right hand side
   double precision       :: ts_coord(3)                   ! stage position 
end type Tytime_stage
!
!AA406 start
!.. defines the structure of the wall elements (BSPH)
type TyParticle_w
   integer(4)             :: cella                         ! cell number 
!AA601   
   integer(4)             :: adjacent_faces(3)             ! index of the adjacent_faces (they are maximum 3)
!AA601 mv   
   integer(4)             :: wet                           ! wet wall particle (distance from fluid particle <= 1.3dx)
   double precision       :: coord(3)                      ! position
   double precision       :: vel(3)                        ! velocity
   double precision       :: dens                          ! density
   double precision       :: pres                          ! pressure
   double precision       :: normal(3)                     ! normal
   double precision       :: weight                        ! area(3D)/length(2D) of the wall element 
   double precision       :: mass                          ! mass of the semi-particles
!AA601 start
   double precision       :: k_d                           ! Coefficient used to estimate the semi-particle volumes
   double precision       :: volume                        ! Semi-particle volume (area in 2D)
!AA601 end
!AA406 end
!
end type TyParticle_w
!
!AA406 end
!
!AA501b start

!.. defines the structure of the body elements
type body_element
   integer(4)             :: npart                                ! number of body particles 
   double precision       :: L_geom(3)                            ! element dimensions (Lx,Ly,Lz): a paralelepiped 
   double precision       :: dx(3)                                ! body particle spacing  
   double precision       :: x_CM(3)                              ! position of the centre of mass
   double precision       :: alfa(3)                              ! rotation angle of the element with respect to the reference 
!                                                                   system
   integer(4)             :: normal_act(6)                        ! boolean value to activate normals if the element side a is 
!                                                                   body side
   double precision       :: mass_deact(6)                        ! (xmin,xmax,ymin,xmax,zmin,zmax) to deactivate particle masses
end type body_element

!.. defines the structure of the bodies
type body
   integer(4)             :: npart                                ! number of body particles 
   double precision       :: mass                                 ! mass
   double precision       :: Ic(3,3)                              ! Moment of inertia
   double precision       :: Ic_inv(3,3)                          ! Inverse of the moment of inertia 
   integer(4)             :: Ic_imposed                           ! Flag to impose Ic in input 
   double precision       :: x_CM(3)                              ! position of the centre of mass
   double precision       :: alfa(3)                              ! rotation angle of the body with respect to the ref. system
   double precision       :: x_rotC(3)                            ! centre of rotation just to configure the initial geometry
   double precision       :: u_CM(3)                              ! velocity of the centre of mass
   double precision       :: omega(3)                             ! angular velocity 
   double precision       :: Force(3)                             ! Force
   double precision       :: Moment(3)                            ! Angular moment
   double precision       :: umax                                 ! Maximum of the absolute value of the particle velocity
   double precision       :: pmax                                 ! Maximum value of pressure
   integer(4)             :: n_elem                               ! number of body particles 
   type(body_element),dimension(:),allocatable :: elem            ! Array of the elements of the body
   integer(4)             :: imposed_kinematics                   ! Flag for imposed kinematics   
   integer(4)             :: n_records                            ! Number of records for imposed body kinematics
   double precision,dimension(:,:),allocatable :: body_kinematics ! Array for the imposed body kinematics (n_records * 7) 
end type body

!
!.. defines the structure of the body particles
type body_particle
   integer(4)             :: body                          ! body number 
   integer(4)             :: cell                          ! cell number 
   double precision       :: rel_pos(3)                    ! relative position (with respect to the body centre of mass)
   double precision       :: pos(3)                        ! position (with respect to the origin of the reference system)
   double precision       :: vel(3)                        ! velocity
   double precision       :: vel_mir(3)                    ! velocity of mirror particles
   double precision       :: acc(3)                        ! acceleration
   double precision       :: normal(3)                     ! normal
   double precision       :: pres                          ! pressure
   double precision       :: mass                          ! mass 
   double precision       :: area                          ! area (length in 2D) for body particles on the body surface
end type body_particle
!AA501b end
!
!.. defines the structure of the zones
!
type TyZone
   character(8)           :: label                         ! nome dell'area iniziale
   character(4)           :: tipo                          ! PERI, SOUR, FLOW, VELO, OPEN, CRIT, LEVE, TAPI, POOL
   character(1)           :: shape                         ! forma dell'area iniziale
   character(1)           :: bend                          ! tipo di colorazione uniforme o a strisce
   character(2)           :: pressure                      ! assegnazione pressione costante o distribuzione idrostatica
   character(3)           :: move                          ! tipo di moto standard o fisso
   character(1)           :: slip                          ! condizione di parete (solo particelle fisse) f=free slip;  n=no slip
   integer(4)             :: ipool                         ! 
   integer(4)             :: npoints                       ! 
   integer(4)             :: icol                          ! colore delle particelle o n° di strisce verticali
   integer(4)             :: Medium                        ! indice mezzo di appartenenza
   integer(4)             :: npointv                       ! numero dei vertici
   integer(4)             :: Indix(2)                      ! 
   integer(4)             :: limit(2)                      ! indici della prima e dell'ultima particella della zona
!AA504 start
   integer(4)             :: IC_source_type                ! IC fluid particle distribution from reservoir vertices and faces (1) or from Cartesian topography (2)
   integer(4)             :: Car_top_zone                  ! zone describing the Cartesian topography, which contains the reservoir (0 if IC_source_type=1)
   integer(4)             :: plan_reservoir_points         ! Number of points describing the reservoir, if (IC_source_type==2)
   integer(4)             :: nag_aux                       ! rough overestimation of the number of fluid particles in the reservoir, if (IC_source_type==2)
   integer(4)             :: ID_first_vertex               ! ID first vertex in case of Cartesian topography
   integer(4)             :: ID_last_vertex                ! ID last vertex in case of Cartesian topography
   integer(4)             :: dam_zone_ID                   ! ID of the dam zone, related to the eventual reservoir (0 if no dam is present), if (IC_source_type==2) 
   integer(4)             :: dam_zone_n_vertices           ! Number of points describing the horizontal projection of the dam zone, if (dam_zone_ID>1) 
   double precision       :: plan_reservoir_pos(4,2)       ! Plane coordinates of the points (3 or 4), which describe the reservoir, if (IC_source_type==2)
   double precision       :: dam_zone_vertices(4,2)        ! Plane coordinates of the points (3 or 4), which describe the dam zone, if (dam_zone_ID>1)
   double precision       :: dx_CartTopog                  ! Spatial resolution of the regular Cartesian Topography 
   double precision       :: H_res                         ! Height of the reservoir free surface
!AA504 end
   double precision       :: pool                          ! 
   double precision       :: coordMM(3,2)                  ! coordinate dei vertici
   double precision       :: vlaw(0:3,MAXPOINTSVLAW)       ! velocita' iniziale
   double precision       :: vel(3)                        ! velocita' iniziale
   double precision       :: trampa                        ! tempo per raggiungere la velocita' della sorgente
   double precision       :: valp                          ! valore pressione iniziale o quota piezometrica
end type TyZone
!
!
!.. defines the structure of the boundary stretch
!
type TyBoundaryStretch
   character(4)           :: tipo                          !["FIXEd", "PERImeter", "POOL", "TAPIs", "SOURce", "LEVEl", "FLOW", "VELOcity", "CRITic", "OPEN"] 
   integer(4)             :: ColorCode                     ! 
   integer(4)             :: numvertices
   integer(4)             :: inivertex
   integer(4)             :: iniside
!  integer(4)              :: numfaces
   integer(4)             :: iniface
   integer(4)             :: medium
   integer(4)             :: zone
   double precision       :: velocity(1:SPACEDIM)          ! velocita' per condizione di TAPI
   double precision       :: NormVelocity                  ! modulo velocita' normale alla faccia
   double precision       :: FlowRate                      ! portata uscente dalla faccia
   double precision       :: trampa                        ! tempo per la rampa della velocita'
   double precision       :: ShearCoeff
   double precision       :: PsiCoeff(1:SPACEDIM)
   double precision       :: FiCoeff(1:SPACEDIM)
end type TyBoundaryStretch
!
!
!.. defines the structure of the medium properties
!
type TyMedium
   character(8)           :: tipo                          ! nome del mezzo ["liquid  ", "gas     ", "general ", "granular", "smagorin"]
   character(8)           :: modelloerosione               ! modello di erosione utilizzato [shields, mohr]
   integer(4)             :: index                         ! numero del mezzo
   integer(4)             :: NIterSol                      ! numero di iterazioni nelle quali il mezzo granular o general è fermo (solido)
   double precision       :: den0                          ! densita' iniziale
!AA504
   double precision       :: den0_s                        ! solid phase reference density (granular flows)
   double precision       :: eps                           ! comprimibilita'
   double precision       :: celerita                      ! celerita' del mezzo
   double precision       :: alfaMon                       ! valore del parametro alfa di Monaghan
   double precision       :: betaMon                       ! valore del parametro beta di Monaghan
   double precision       :: visc                          ! valore della viscosita' cinematica 
   double precision       :: coes                          ! valore della coesione 
   double precision       :: numx                          ! valore della viscosita' cinematica massima
   double precision       :: mumx                          ! valore della viscosita' dinamica massima  (utilizzata per soglia movimento granulare)
   double precision       :: taucri                        ! valore dello sforzo critico
   double precision       :: cuin                          ! esponente della curva reologica
   double precision       :: phi                           ! angolo di attrito pe il materiale granulare
   double precision       :: cons                          ! consistenza   
   double precision       :: Cs                            ! Costante modello Smagorinsky    
   double precision       :: RoughCoef                     ! Coefficiente di rugosita' del mezzo
   double precision       :: D50                           ! Diametro mediano del mezzo granulare
   double precision       :: SettlingCoef                  ! Coefficiente per la velocita' di sedimentazione del mezzo (granulare)
   !!!! modello bifluido
   double precision       :: codif                         ! valore del coefficiente di diffusione
   !!!! modello bifluido
   !! modello esplosioni
   double precision       :: Gamma                         ! state equation constant for explosion
   double precision       :: InitialIntEn                  ! Initial Specific Internal Energy 
   !! modello esplosioni
!AA504 start
   double precision       :: gran_vol_frac_max             ! Maximum volume fraction of the solid granular phase within an SPH granular particle (granular flows)
   double precision       :: d_90                          ! 90-th percentile diameter of the granular size distribution (granular flows)
!AA504 end
end type TyMedium
!
!
!.. defines the structure of the boundary side
!
type TyBoundarySide
   character(4)           :: tipo                           ! ["FIXEd", "PERImeter", "TAPIs", "SOURce", "LEVEl", "FLOW", "VELOcity", "CRITic", "OPEN"] 
   integer(4)             :: stretch
   integer(4)             :: previous_side
   integer(4)             :: vertex(1:SPACEDIM-1)
   integer(4)             :: CloseParticles                 ! numero di particelle vicine al contorno
   double precision       :: length
   double precision       :: CloseParticles_maxQuota        ! massima quota tra le particelle vicine al contorno
!AA601 comm   
   double precision       :: T(1:SPACEDIM, 1:SPACEDIM)      ! Direction cosines of the local reference system of the boundary (the last one is the normal)
   double precision       :: R(1:SPACEDIM, 1:SPACEDIM)      
   double precision       :: RN(1:SPACEDIM, 1:SPACEDIM)
   double precision       :: angle
   double precision       :: velocity(1:SPACEDIM)           ! velocita'
end type TyBoundarySide
!
!
!.. defines the structure of the nodes
!
type TyNode
   integer(4)             :: name
   double precision       :: GX(1:SPACEDIM)
   double precision       :: LX(1:SPACEDIM)
end type TyNode
!
!
!.. defines the structure of the faces
!
type TyBoundaryFace
   type(TyNode)           :: Node(1:MAXFACENODES)
   integer(4)             :: nodes
   integer(4)             :: stretch
   integer(4)             :: CloseParticles                ! numero di particelle vicine al contorno
   double precision       :: area
   double precision       :: CloseParticles_maxQuota       ! massima quota tra le particelle vicine al contorno
   double precision       :: T   (1:SPACEDIM,1:SPACEDIM)   ! Direction cosines of the local reference system of the boundary (the last one is the normal)
   double precision       :: RPsi(1:SPACEDIM,1:SPACEDIM)
   double precision       :: RFi (1:SPACEDIM,1:SPACEDIM)
   double precision       :: velocity(1:SPACEDIM)          ! velocita'
end type TyBoundaryFace
!
!
!.. defines the structure of the close boundary data table
!
type TyBoundaryData
   integer(4)             :: CloBoNum                      ! numero facce di contorno che interessano una particella
   double precision       :: LocXYZ(1:SPACEDIM)            ! normali alla faccia in coordinate locali
   double precision       :: BoundaryIntegral(1:8)         ! tabella integrali 
   double precision       :: IntGiWrRdV(1:SPACEDIM,1:SPACEDIM) ! tabella integrali 
end type TyBoundaryData
!
!
!.. defines the structure of the  
!
type TyBoundaryConvexEdge
   type(TyNode)           :: Node(1:2)                     !
   integer(4)             :: face(1:2)                     ! 
   double precision       :: length                        ! 
   double precision       :: component(1:SPACEDIM)         ! 
end type TyBoundaryConvexEdge
!
!
!.. defines the structure of the control points
!
type TyCtlPoint
   integer(4)             :: cella                         ! 
   double precision       :: coord(3)                      ! coordinate
   double precision       :: vel(3)                        ! velocita' interpolata nel punto
   double precision       :: pres                          ! pressione interpolata nel punto
   double precision       :: dens                          ! densita' interpolata nel punto
   double precision       :: uni                           ! unita' interpolata nel punto
   double precision       :: dist                          ! 
end type TyCtlPoint
!
!
!.. defines the structure of the sections
!
type TySection
   character(8)           :: Label                         ! 
   character(2)           :: Tipo                          ! x,y,x,g
   integer(4)             :: Icont(2)                      ! puntatori ai punti di controllo generati
   integer(4)             :: ColorCode                     !    ???
   double precision       :: Constant(1:SPACEDIM)          ! 
   double precision       :: XYZRange(1:SPACEDIM,2)        ! 
   double precision       :: TGLsection(1:SPACEDIM,1:SPACEDIM)  !
end type TySection
!
!
!.. defines the structure of the control lines
!
type TyCtlLine
   character(8)           :: label                         ! nome della linea
   integer(4)             :: icont(2)                      ! puntatore iniziale e finale nell'insieme dei punti
end type TyCtlLine
!

!AA504 start
 type tyQ_section_array                                             ! Section ID depends on the input order
    integer(4)                                 :: n_vertices        ! number of vertices describing a monitoring section for the flow rate (3 or 4)
    double precision,dimension(:),allocatable  :: flow_rate         ! section flow rate per fluid type and global: flow_rate(n_fluid_types+1) 
    double precision                           :: area              ! section area 
    double precision                           :: normal(3)         ! section normal 
    double precision                           :: vertex(4,3)       ! section vertices (in case of 3 vertices do not mind about the fourth point)
    double precision                           :: vertex_loc(4,3)   ! section vertices (local coordinates: the origin is vertex1, the third local axis is the normal)
    double precision                           :: loc_axis(3,3)     ! local axis (axis 1: vertex2-vertex1, axis 2: vector product of local axis 1 and normal; axis 3: normal)
 end type
 
 type TyQ_section
    integer(4)                                         :: n_sect            ! number of monitoring sections for the flow rate 
    integer(4)                                         :: n_fluid_types     ! number of fluid types (the first ID fluid types are selected)
    integer(4)                                         :: it_out_last       ! auxiliary variable to print flow rates 
    double precision                                   :: dt_out            ! writing time step for section flow rates
    type(tyQ_section_array), dimension(:), allocatable :: section           ! type section to monitor flow rates
 end type  

 type TyGranular_flows_options 
    integer(4)                                         :: ID_erosion_criterion    ! erosion criterion ID (1:Shields-Seminara,2:Shields,3:Mohr-Coulomb)
    integer(4)                                         :: ID_main_fluid           ! ID of the main flow 
    integer(4)                                         :: ID_granular             ! ID of the granular medium
    integer(4)                                         :: monitoring_lines        ! number of monitoring lines aligned with x- or y-axis
    integer(4)                                         :: it_out_last             ! auxiliary variable for printing
    integer(4)                                         :: n_max_iterations        ! maximum number of iterations (iterative process, which estimates u_star)
    integer(4)                                         :: erosion_flag            ! erosion_flag: 0(activated far from fronts); 1(not activated), 2(activated everywhere)
    integer(4)                                         :: viscosity_blt_formula   ! viscosity_blt_formula: formula for viscosity in the bed load transport region 
                                                                                  ! (1:Chauchat-Medale 2010 CMAME; 2:Chezy-like; 3:diluted viscosity)
    integer(4)                                         :: deposition_at_frontiers ! forced deposition at frontiers (erosion criterion=2): yes(1) or no(2) 
    integer(4)                                         :: Gamma_slope_flag        ! flag to activate (or not) Gamma_slope (effects only when ID_erosion criterion=1) 
    double precision                                   :: conv_crit_erosion       ! convergence criterion for erosion  
    double precision                                   :: velocity_fixed_bed      ! (optional) velocity_fixed_bed: velocity threshold (input) to speed-up a simulation 
    double precision                                   :: dt_out                  ! dt_output for the interface monitoring lines
    double precision                                   :: Chezy_friction_coeff    ! Chezy's friction coefficient (in case of Chezy-like viscosity formula)
    double precision                                   :: x_min_dt                ! x_min to involve SPH granular particles in dt estimation
    double precision                                   :: x_max_dt                ! x_max to involve SPH granular particles in dt estimation
    double precision                                   :: y_min_dt                ! y_min to involve SPH granular particles in dt estimation
    double precision                                   :: y_max_dt                ! y_max to involve SPH granular particles in dt estimation
    double precision                                   :: z_min_dt                ! z_min to involve SPH granular particles in dt estimation
    double precision                                   :: z_max_dt                ! z_max to involve SPH granular particles in dt estimation    
    double precision,dimension(:,:),allocatable        :: lines                   ! x and/or y coordinates defining the monitoring line/point
 end type 
!AA504 end

!AA601 start
type face_der_type                                                                 ! Used by DB-SPH 
   integer(4),dimension(4)                             :: vert_list                ! List of face/segment vertices (triangles in 3D)  
   double precision                                    :: area                     ! Face area/length in 3D/2D
   double precision,dimension(3)                       :: normal                   ! Face normal 
end type

type vertex_der_type
   double precision,dimension(3)                       :: pos                      ! Vertex position                                         
end type

type DBSPH_surf_mesh_der_type 
   type(vertex_der_type),allocatable,dimension(:)      :: vertices                 ! List of vertices of the surface mesh 
   type(face_der_type),allocatable,dimension(:)        :: faces                    ! List of faces (both 3D and 2D)
end type 

type DBSPH_der_type 
!AA601 DBSPH%n_w replaces nag_w everywhere    
   logical                                             :: MUSCL_boundary_flag      ! Flag to activate (or not) boundary terms for MUSCL reconstruction
   logical                                             :: in_built_monitors        ! Flag to activate (or not) in-built motion of control lines 
   integer(4)                                          :: n_w                      ! Number of surface elements
   integer(4)                                          :: n_monitor_points         ! Number of monitoring points
   integer(4)                                          :: n_monitor_regions        ! Number of monitoring regions (0 or 1) to estimate the Force along x-direction 
   integer(4)                                          :: n_kinematics_records     ! Number of records for imposed kinematics
   integer(4)                                          :: n_inlet                  ! Number of inlet sections (to impose DBSPH inlet BC)
   integer(4)                                          :: n_outlet                 ! Number of outlet sections (to impose DBSPH outlet BC)
   double precision                                    :: dx_dxw                   ! Ratio between the fluid and the semi-particle sizes
   double precision                                    :: k_w                      ! Coefficient to compute semi-particle volumes (ref. DBSPH equations)
   double precision                                    :: monitor_region(6)        ! (xmin,xmax,ymin,xmax,zmin,zmax) to detect the monitoring region
   type(DBSPH_surf_mesh_der_type)                      :: surf_mesh                ! DBSPH surf mesh: vertices and faces 
   integer(4),allocatable,dimension(:)                 :: monitor_IDs              ! IDs of the monitoring points    
   double precision,dimension(:,:),allocatable         :: kinematics               ! Array for the imposed kinematics of the DBSPH surface elements (n_records,4)

   double precision,dimension(:,:),allocatable         :: inlet_sections           ! Array for the inlet sections useful to DBSPH inlet BC (n_records,10): 
!                                                                                    positions, normal, velocity, length
   double precision,dimension(:,:),allocatable         :: outlet_sections          ! Array for the outlet sections useful to DBSPH outlet BC (n_records,7): 
!                                                                                    positions, normal, length, pressure
end type 
!AA601 end

!------------------------------------------------------------------
!.. variables for domain, grid and set to zero the current particle array
type (TyGlobal)     :: Domain
type (TyGriglia)    :: Grid
type (TyParticle)   :: PgZero
!AA402
type (Tytime_stage) :: ts_pgZero
!AA504
type (TyQ_section)  :: Q_sections
type (TyGranular_flows_options) Granular_flows_options
!AA601
type(DBSPH_der_type) :: DBSPH
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. variabile per convex edges
!AA504 rm
!type(TyBoundaryConvexEdge),dimension(:),allocatable :: BoundaryConvexEdge
!------------------------------------------------------------------

end module AdM_User_Type
