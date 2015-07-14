!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ALLOC_Module
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Add the message unit on the screen
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Amicarelli/Agate  05/07/2011     Time stage parameters
! 04  Amicarelli/Agate  30Nov11        BSPH: wall element parameters
!
!AA501b
! 05  Amicarelli-Agate  13Nov12        Body dynamics
!AA504
! 06  Amicarelli        8Apr14         (v5.04) Variables for some memory management from input file
!
!************************************************************************************
! Module purpose : allocatable arrays declaration
!
! Calling routine: all
!
! Called routines: none
!
!************************************************************************************
!
module ALLOC_Module

use AdM_User_Type

!------------------------------------------------------------------
! allocatable arrays
!------------------------------------------------------------------
!type (TyVertex),          dimension(:),   allocatable :: Vertice

!AA504 sub comm
double precision,         dimension(:,:), allocatable :: Vertice                  ! Vertice (3,n_vertices): array of the geometrical vertices

type (TyBoundaryStretch), dimension(:),   allocatable :: Tratto
type (TyZone),            dimension(:),   allocatable :: Partz
type (TyMedium),          dimension(:),   allocatable :: Med
integer(4),               dimension(:),   allocatable :: BoundaryVertex
type (TyBoundarySide),    dimension(:),   allocatable :: BoundarySide
type (TyBoundaryFace),    dimension(:),   allocatable :: BoundaryFace
type (TyCtlPoint),        dimension(:),   allocatable :: Control_Points
type (TyCtlLine),         dimension(:),   allocatable :: Control_Lines
type (TySection),         dimension(:),   allocatable :: Control_Sections
type (TyCtlPoint),        dimension(:),   allocatable :: Section_Points
type (TyParticle),        dimension(:),   allocatable :: Pg
!AA402
type (Tytime_stage),      dimension(:),   allocatable :: ts0_pg
!AA504
type(TyBoundaryConvexEdge),dimension(:),  allocatable :: BoundaryConvexEdge
!
!AA504sub comment
!AA406
! Icont(grid%nmax+1) contains the number of the first particle in each cell and the total number of particles.
!    The particle are here ordered according to the cell order. The element values of this array monotonically increase. 
! NPartOrd(PARTICLEBUFFER) contains the particle IDs. The elements are ordered as for Icont.
integer(4),               dimension(:),   allocatable :: Icont, NPartOrd, GCBFVector
!
!AA406 start
! Array of the wall elements (BSPH)
type (TyParticle_w),      dimension(:),   allocatable :: Pg_w
! Arrays for ordering wall elements, according to the underliying mesh
integer(4),               dimension(:),   allocatable :: Icont_w, NPartOrd_w
!AA406 end

!AA501b start
! Array of the transported rigid solid bodies
type (body),              dimension(:),   allocatable :: body_arr
! Array of the body particles
type (body_particle),     dimension(:),   allocatable :: bp_arr
! Arrays for ordering body particles, according to the underliying mesh
integer(4),               dimension(:),   allocatable :: Icont_bp, NPartOrd_bp
!AA501b end

integer(4),               dimension(:,:), allocatable :: GCBFPointers
integer(4),               dimension(:),   allocatable :: BFaceList
!------------------------------------------------------------------
!
!AA406 sub
! nPartIntorno(PARTICLEBUFFER): array of the number of the neighbouring particles
!
integer(4),               dimension(:),   allocatable :: nPartIntorno
!
!AA406
! PartIntorno(NMAXPARTJ*PARTICLEBUFFER): array of the indeces of the neighbouring particles
!
integer(4),               dimension(:),   allocatable :: PartIntorno
!
!AA406 comm: PartKernel(4,NMAXPARTJ*PARTICLEBUFFER)
!            PartKernel(1,b): -|gradW_0b|/|r_0b|, gradW: kernel gradient (cubic spline); 
!                             Then gradW = |gradW_0b| * (gradW_0b/|gradW_0b|) = - PartKernel*rag, 
!                             rag=-(x_b-x_0) è orientato come gradW
!AA504 sub comment
!            PartKernel(2,b): -|gradW_0b|/[|r_0b|(|r_0b|^2+eps^2)], cubic spline kernel
!            PartKernel(3,b): -|gradW_0b|/|r_0b|, Gallati anti-cluster kernel, used for pressure terms (by "semi" approach)
!            PartKernel(4,b): W_0b: module of the kernel (cubic spline), used for interpolations and BSPH
! gradW vector is equal to -PartKernel(1 or 3,b)*rag_0b
!
double precision,         dimension(:,:), allocatable :: PartKernel
!
!AA406 comm: rag(3,NMAXPARTJ*PARTICLEBUFFER): 3D vector list of -r_0b=x_0-x_b,
!            r_0b (vector distance between the computational particle and the neighbour)
!
double precision,         dimension(:,:), allocatable :: rag  
!
!AA406 start
! nPartIntorno_fw(PARTICLEBUFFER): array of the number of the neighbouring wall particles 
integer(4),               dimension(:),   allocatable :: nPartIntorno_fw
! PartIntorno_fw(NMAXPARTJ*PARTICLEBUFFER): array of the indeces of the neighbouring wall particles 
integer(4),               dimension(:),   allocatable :: PartIntorno_fw   
! Kernel function neighbouring array (wall neighbours), kernel_fw(NMAXPARTJ*PARTICLEBUFFER)
!AA406test
double precision,         dimension(:,:),   allocatable :: kernel_fw
! relative distances from wall particles: -r_0a, rag_fw(components,NMAXPARTJ*PARTICLEBUFFER)
double precision,         dimension(:,:), allocatable :: rag_fw  
!AA406 end

!AA501b start
! neighbouring arrays for body dynamics
! neighbouring arrays of the body particles (body particle - fluid particle interactions; 
!                                            fluid particles are here treated as neighbours)
! nPartIntorno_bp_f(n_body_particles): array of the number of the neighbouring fluid particles to each body particle
integer(4),               dimension(:),   allocatable :: nPartIntorno_bp_f
! PartIntorno_bp_f(NMAXPARTJ*n_body_particles): array of the indeces of the neighbouring fluid particles to each body particle
integer(4),               dimension(:),   allocatable :: PartIntorno_bp_f   
! Kernel derivative neighbouring array (fluid neighbours), KerDer_bp_f_cub_spl(n_body_particles*NMAXPARTJ), cubic spline;
! KerDer_bp_f_cub_spl  = -|gradW_bp_f|/|r_bp_f|
double precision,         dimension(:), allocatable :: KerDer_bp_f_cub_spl 
! Kernel derivative neighbouring array (fluid neighbours), KerDer_bp_f_Gal(n_body_particles*NMAXPARTJ), Gallati's kernel;
! KerDer_bp_f_Gal  = -|gradW_bp_f|/|r_bp_f|
double precision,         dimension(:), allocatable :: KerDer_bp_f_Gal 
! relative distances from fluid particles: -r_bp_f; rag_bp_f(3,NMAXPARTJ*n_body_particles)
double precision,         dimension(:,:), allocatable :: rag_bp_f  
! neighbouring arrays for inter-body impacts (body_A particle - body_B particle interactions)
! list of surface body particles surf_body_part(n_surf_body_part)
integer(4),               dimension(:),   allocatable :: surf_body_part
! nPartIntorno_bp_bp(n_surf_body_part): array of the number of the neighbouring body particles, belonging to another body
integer(4),               dimension(:),   allocatable :: nPartIntorno_bp_bp
! PartIntorno_bp_bp(n_surf_body_part*NMAXPARTJ): array of the indeces of the neighbouring body particles (of another body) 
integer(4),               dimension(:),   allocatable :: PartIntorno_bp_bp   
! relative distances from body particles (belonging to another body): -r_bp_bp; rag_bp_bp(3,NMAXPARTJ*n_surf_body_part)
double precision,         dimension(:,:), allocatable :: rag_bp_bp  
! array of velocity impacts for body dynamics impact_vel(n_surf_body_part x (n_bodies+n_boundaries))
double precision,         dimension(:,:), allocatable :: impact_vel  
!AA501b end

!------------------------------------------------------------------
!------------------------------------------------------------------
!.. arrays to calc the integrals table
integer(4),               dimension(:,:), allocatable :: BoundaryDataPointer
type (TyBoundaryData),    dimension(:),   allocatable :: BoundaryDataTab
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. arrays to count in and outgoing particles for every medium
integer(4),               dimension(:),   allocatable :: OpCount,SpCount,EpCount,EpOrdGrid
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. array to count fluid particles (not 'sol' in case of medium 'granular')
integer(4),               dimension(:),   allocatable :: Array_Flu
!integer(4),               dimension(:),   allocatable :: Array_Sol
!------------------------------------------------------------------
!------------------------------------------------------------------
!.. arrays 2D to calc the free surface in case of erosion model
integer(4),               dimension(:,:,:),allocatable :: ind_interfaces    ! store pl_imed, intliq_id, intsol_id, of every colomn of cells
!AA504 rm
!------------------------------------------------------------------

!AA504
real(kind=kind(1.d0)),dimension(:),allocatable :: Z_fluid_max,q_max  ! (only in 3D) The 2D arrays of the maximum values of the fluid particle height (at the 
                                                               ! nodes of the grid vertical columns) and the specific flow rate 

end module
