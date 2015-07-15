!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : GLOBAL_Module
!
! Last updating : May 02, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2011           varie
! 04  Amicarelli        23/11/2011     multiple inlet
! 05  Amicarelli/Agate  30Nov11        BSPH: wall elements and parameters
! 06  Agate             07/05/2012     Limited and licensed version
! 07  Agate             15/05/12       err file for erosion model
! 08  Agate             03/07/12       Module exclusion (Diffusion, Esplosion, Temporal schemes)
! 09  Amicarelli-Agate  13nov12        (AA501b) Body dynamics       
! 10  Agate             27/05/14       Add more modules check in license
!
!************************************************************************************
! Module purpose : Global variables declaration
!
! Calling routine: all
!
! Called routines: none
!
!************************************************************************************
!
module GLOBAL_Module
!
!------------------------------------------------------------------
character(len=8),  parameter :: acode     = "Sphera  "
!
character(len=8),  parameter :: version   = "6.0.0   "
!
character(len=14), parameter :: date_vers = "May       2014"
character(len=8)             :: dt_opt    = "original"   !"euristic"
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. Sphera Limited version: the limitation to PARTICLE_NUMBER_LIMIT number of particles
!..  =0 no limits
!..  >0 max number of particles
integer(4), public :: PARTICLE_NUMBER_LIMIT  ! codice L numeroparticelle
!------------------------------------------------------------------
!.. Modules excluded in Sphera
!..  = .true.   the module is available in Sphera
!..  = .false.  the module is not available in Sphera
logical, public    :: Erosion_Module_Shields_Mohr        ! codice  1     1 -    1
logical, public    :: Diffusion_Module                   ! codice  2     2 -    3
logical, public    :: Explosion_Module                   ! codice  3     4 -    7
logical, public    :: TemporalScheme_Module              ! codice  4     8 -   15
logical, public    :: BodyDynamics_Module                ! codice  5    16 -   31
logical, public    :: DBSPH_Module                       ! codice  6    32 -   63
logical, public    :: MultiFluid_Module                  ! codice  7    64 -  127
logical, public    :: MoreFluids_Module                  ! codice  8   128 -  255
logical, public    :: Granular_flux                      ! codice  9   256 -  511
logical, public    :: Erosion_Module_Shields_Seminara    ! codice 10   512 - 1023
!------------------------------------------------------------------
!Proposta giu2013 moduli Sphera ("no" da togliere, "si" da mettere):
!1-    erosion
!6-    DBSPH (AA)
!4-    Time integration (AA)
!5-    Body (AA)
!no-    threshold: n particles
!no7-    multi-fluid (no validation)
!9-    viscous flows (no validation)
!3-    explosion (no validation)
!2-    diffusion (no validation)
!si-    SASPH
!si-    Visco-elastic pseudo-Newtonian
!------------------------------------------------------------------
integer(4), parameter :: MAXTIT = 10
!
character(80), dimension(MAXTIT) :: title
!------------------------------------------------------------------
!------------------------------------------------------------------
!AA504 sub
integer(4), public, parameter :: lencard= 200
!
character(255) :: nomecaso, nomecas2
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. pointer for coordinates location 2d or 3d
integer(4),    dimension(0:3,2) :: icoordp    ! = (/ 0, 1,3,0,  0, 1,2,3 /)  inizializzato nel main per compatibilita con xlf90
character( 1), dimension(0:3)   :: xyzlabel   = (/ "T", "X", "Y", "Z" /)  
character( 4), dimension(3)     :: ncordlabel = (/ "    ", "(2D)", "(3D)" /)  
!
!.. maximum number of points for definition of GENERIC area
integer(4),       public, parameter :: MAXPOINTSZONE = 20
!maximum number of points for definition of velocity LAW
integer(4),       public, parameter :: MAXPOINTSVLAW  = 50  !10
!
!.. global constants for arrays dimension
integer(4),       public, parameter :: SPACEDIM = 3
integer(4),       public, parameter :: PLANEDIM = 2
!AA504 sub
integer(4),       public, parameter :: MAXCLOSEBOUNDSIDES = 2 !v503: 2   !5
!AA504 rm line
integer(4),       public, parameter :: MAXOPENSIDES = 10
integer(4),       public, parameter :: MAXOPENFACES = 50
integer(4),       public, parameter :: MAXPARTLINE = 2000
integer(4),       public, parameter :: LIMCLOSEBOUNDSIDES = 2
integer(4),       public, parameter :: MAXFACENODES = 4
!AA405 sub
!AA504 rm line
integer(4),       public, parameter :: NUMCOLS_BIT = 5
integer(4),       public, parameter :: NUMROWS_BIT = 22
integer(4),       public, parameter :: INT_KERNELTABLE = 40
!AA405 sub
integer(4),       public, parameter :: INIPARTICLEBUFFER = 225000 !100000 (altrimenti non abbastanza particelle per getti)
!
!integer(4),       public, parameter :: FI_STEPS = 5
!integer(4),       public, parameter :: TETA_STEPS = 20
!integer(4),       public, parameter :: FI_STEPS = 10
integer(4),       public, parameter :: FI_STEPS = 8
integer(4),       public, parameter :: TETA_STEPS = 4 * FI_STEPS
!
!.. global constants
integer(4),       public, parameter :: menouno = -1
double precision, public, parameter :: max_positive_number = 3.4d+38
double precision, public, parameter :: max_negative_number = -3.4d+38
double precision, public, parameter :: PIGRECO = 3.14159265358979              ! pi-greco as constant
!double precision, public, parameter :: MAXUNI = 0.95d0
double precision, public, parameter :: KERNELCONST2D = 0.454728408833987d0     ! =10/(7*PIGRECO)
double precision, public, parameter :: KERNELCONST3D = 1.0d0 / PIGRECO         ! 1/pigreco as constant
double precision, public, parameter :: FI_INTERVAL = PIGRECO/2.0d0             ! =PIGRECO/2 
double precision, public, parameter :: TETA_INTERVAL = PIGRECO + PIGRECO       ! =2*PIGRECO
double precision, public, parameter :: zero    = 0.0d0                         ! zero as constant
double precision, public, parameter :: one     = 1.0d0                         ! one as constant
double precision, public, parameter :: two     = 2.0d0                         ! two as constant
double precision, public, parameter :: four    = 4.0d0                         ! four as constant
double precision, public, parameter :: half    = 0.5d0                         ! 1/2 as constant
double precision, public, parameter :: quarter = 0.25d0                        ! 1/4 as constant
double precision, public, parameter :: xyz_tolerance= 1.0d-4                   ! tolerance for coord checking
double precision, public, parameter :: const_m_9999 = -9999.0d0                ! -9999. as constant
double precision, public, parameter :: sqrttwo = 1.4142135623731d0             ! square root of two as constant
double precision, public, parameter :: arrotondamento = 1.0d-5                 ! tolerance
double precision, public, parameter :: azzeramento = 1.0d8                     ! tolerance
double precision, public, parameter :: GI = 9.80665                            ! gravitational constant
double precision, public, parameter :: vKconst = 0.41                          ! von Karman constant
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. global variables 
!AA601
logical          :: dt_alfa_Mon
integer(4)       :: Ncord, NMedium, NPartZone
integer(4)       :: NPointst,NPoints,NPointsl,NPointse, NLines, NSections, GCBFVecDim
integer(4)       :: NumVertici, NumFacce, NumTratti, NumBVertices, NumBSides
!
!AA601 sub
integer(4)       :: Nag,nagpg,PARTICLEBUFFER,nag_reservoir_CartTopog
!
integer(4)       :: it_start, it_corrente, it_eff
integer(4)       :: indarrayFlu, indarraySol
double precision :: tempo, dt, pesodt, dt_average, DTminBER
!
!.. variables to check the cpu time in parallel version
!!!double precision :: cpu_prima,cpu_dopo, cpu_loop1a,cpu_loop1b,cpu_loop2a,cpu_loop2b,cpu_loop3a,cpu_loop4a,cpu_loop5a,cpu_loop6a
!
!.. variabiles for dimensioning the array that store the gradients of the particles
!..  surrounding the current one and the array that store the close boundaries and the integrals
!.. for the current particle in every loop
!
!AA503sub and rm
!double precision, public, parameter :: COEFNMAXPARTJ = 1.5d0 !Sphera 2.0d0 !1.0d0  !1.5d0  !explosion 4.0d0  !original 1.5d0

!AA504rm line
integer(4)     :: NMAXPARTJ            ! max number of particles surrounding the current one
integer(4)     :: MaxNcbs              ! max number of close boundary sides for the current particle
integer(4)     :: MaxNcbf              ! max number of close boundary faces for the current particle
!
character(255) :: LicenseFile          ! file name for license 
character(255) :: nomefilekill         ! file name for kill the execution
!AA504
character(10)  :: exetype              ! type of machine: "windows" or "linux"
logical        :: kill_flag            ! flag for kill the execution
character(255) :: nomefileerr          ! file name for err file in erosion model
logical        :: err_flag             ! flag for err file existence in erosion model
logical        :: Restart              ! flag if this run is a restart
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. general purpose variables
! 
  double precision :: eta              ! equivalent to 0.001 * Domain%h
  double precision :: eta2             ! equivalent to 0.01 * Domain%h * Domain%h
  double precision :: doubleh          ! equivalent to 2.*Domain%h
  double precision :: squareh          ! equivalent to Domain%h*Domain%h
  double precision :: doublesquareh    ! equivalent to doubleh * doubleh
  double precision :: cubich           ! equivalent to Domain%h*Domain%h*Domain%h
  double precision :: unosuh           ! equivalent to 1./Domain%h
  double precision :: unosusquareh     ! equivalent to 1./(Domain%h*Domain%h)
!
  logical          :: current_version  ! 
  logical          :: diffusione       ! flag if the diffusion model (bifluid) is invoked
  logical          :: esplosione       ! flag if the explosion model is invoked
  logical          :: erosione         ! flag if the granular for erosion model is invoked
  character(8)     :: modelloerosione  ! type of erosion model is invoked (shields, mohr)

!AA501b start
  integer(4)       :: n_body_part      ! total number of body particles
  integer(4)       :: n_surf_body_part ! total number of body particles belonging to the body surfaces
  integer(4)       :: n_bodies         ! total number of bodies
  integer(4)       :: imping_body_grav ! flag to consider (1) or not (0) gravity components perpendicular &
                                       ! to the normal vectors of the impingement surfaces
  double precision :: dx_dxbodies      ! ratio between fluid particle and body particle size   
!AA501b end
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. global variables for source particles (eliminazione problema opzione compilatore AUTOMATIC)
integer(4)       :: NumOpenSides, SourceSide, mat, irz, izone
!
!AA406 sub
integer(4)       :: SourceFace, NumOpenFaces,pinttimeratio, itime_jet
!
!AA405 rm
!double precision :: RowVelocity, RowPeriod, yfila, ParticleVolume, zfila
!
!AA405 start
double precision :: RowPeriod, yfila, ParticleVolume, zfila
double precision, dimension(1:MAXOPENSIDES)      :: RowVelocity
integer(4), dimension(1:MAXOPENSIDES)            :: NumPartperLine, NumPartFace
!AA405 end
!
integer(4),      dimension(1:MAXOPENSIDES)       :: Openside
integer(4),      dimension(1:MAXOPENFACES)       :: OpenFace
!
!AA405 sub
double precision,dimension(1:MAXOPENSIDES, 1:MAXPARTLINE, 1:SPACEDIM) :: PartLine
!
double precision,dimension(1:SPACEDIM)           :: P
double precision,dimension(1:SPACEDIM)           :: Q
double precision,dimension(1:SPACEDIM)           :: nn
!------------------------------------------------------------------
!.. global variables for integrals calculation
integer(4),parameter                             :: BITrows = FI_STEPS * TETA_STEPS
!integer(4),parameter                             :: BITcols = 4
integer(4),parameter                             :: BITcols = 7
integer(4),parameter                             :: ktrows = INT_KERNELTABLE
integer(4),parameter                             :: ktcols = 4
double precision                                 :: ktdelta
double precision,dimension(0:ktrows, 0:ktcols)   :: kerneltab
double precision,dimension(1:BITrows, 1:BITcols) :: BoundIntegralTab
!------------------------------------------------------------------
!.. global variables for vtkconverter
integer(4), public, parameter           :: maxnumblock = 9999   ! maximum number of blocks (memorizations)
integer(4)                              :: nblocchi,block       ! number of current block 
double precision                        :: freq_time, val_time  ! frequency time to write the blocks
logical                                 :: vtkconv              ! flag if the vtkconverter is invoked
integer(4),      dimension(maxnumblock) :: blocchi              ! array to store the number of blocks
double precision,dimension(maxnumblock) :: Time_Block           ! array to store the time that the block is stored
!------------------------------------------------------------------
!.. global variables for diffusion model 
  double precision, public, parameter   :: VFmn = 0.1d0         ! frazione di volume minima
  double precision, public, parameter   :: VFmx = 0.9d0         ! frazione di volume massima
!------------------------------------------------------------------
!.. global variables for find number of edges
  integer(4)                            :: NumBEdges            ! number of Edges (numero di spigoli convessi)
!------------------------------------------------------------------
!.. global variables for indices of celle that must be considered around the current one in sub CalcVarLength
  integer(4), dimension(14,3)           :: indicecelle          ! sequenza scansione celle con un nuovo sistema in CalcVarLength
!------------------------------------------------------------------
!
end module
