!cfile sphera.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : sphera
!
! Last updating : April 08, 2014
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2011           varie
! 04  Amicarelli/Agate  30/06/2011     Time stage parameters
! 05  Amicarelli        23/11/2011     multiple inlet
! 06  Agate             07/05/2012     Limited and licensed version
! 07  Agate             15/05/2012     err file for erosion model
! 08  Amicarelli-Agate  13nov12        Body dynamics
!AA504 start
! 09  Amicarelli        08Apr14        (v5.04) 
!                                      3D 2-phase model for granular flows: 
!                                         *SPH model based on SPH approximation of Chauchat-Médale 2D FV model (2010, CMAME)
!                                         *Other modifications:
!                                            Correction in estimating I2(2e_ij) -second invariant of the strain rate tensor- in 3D
!                                            Correction of the velocity gradients in 3D (now without renpormalization)
!                                            I2(2e_ij)  computed only when necessary
!                                            Options for viscosity formulation: Chauchat-Médale mixture viscosity (no tuning parameters), Chezy-like viscosity 
!                                                                               (with a tuning parameter): diluted viscosity. 
!                                            Optional and faster solution with the fixed bed detected according to a velocity threshold
!                                            No-slip conditions in 3D are disabled (lack of validation for the boundary term depending on molecular viscosity -mom.equation-)
!                                            Neglection of the moleclar viscosity term for Stokes compressible fluids in momentum equation 
!                                            Monaghan's viscosity always active (also for separating elements)
!                                      Erosion criterion:
!                                         *Shields - van Rijn - Seminara et al. 3D erosion criterion
!                                         *Granular mixture can erode the fixed bed (physical depth of the bed load transport region) 
!                                         *Other modifications:
!                                            k(roughness)=3*d_90 for 3D erosion criterion
!                                            Correct re-initialization of SPH granular particles
!                                            Main bed slope angle computed according to the interface normal in 2D
!                                            Shields threshold for low Re* according to data of Dey (1999, AMM)
!                                            Convergence criterion and maximum iteration numer now in input file
!                                            Corrections in the definition of the pure fluid - mixture interface 
!                                      *Numerics:
!                                         Some memory optimization (some memory management according to input data)
!                                         Optimizing CPU time for initiialization of water reservoirs and intersections between the underlying grid and boundary faces 
!                                         Alternative tretment of reservoirs, provided as PV topographic file + pre-processing in Sphera
!                                         Post-processing for: hydrographs, water depth, specific flow rate, time evolution of free surface and mixture-pure fluid interface 
!                                         Corrected parallelization (no dependency on threat order) in defining the interfaces (granular flows and erosion criterion)
!                                         Time step depends on mixture viscosity and Monaghan viscosity too
!                                         Time only depends on SPH granular particles free to move 
!AA504 end
!AA601 start
! 10  Amicarelli        26Jan15        DBSPH-input(AA601): main features (Amicarelli et al. 2013, IJNME)
!                                         Subroutines to import and process surface meshes from SHM for the DBSPH scheme (SHM-OpenFoam - Sphera coupling).
!                                         New formula for the DBSPH shape coefficient of the surface elements (semi-particle volumes), 
!                                            based on the angle formed between adjacent surface elements. 
!                                         Initialization of the integral Shepard coefficient (for fluid particles) using fictitious and temporary reservoirs.
!                                         DB-SPH hard-coding integrated in the main code.
!                                      New subroutines:    
!                                         Import_ply_surface_meshes, ReadDBSPH, DBSPH_IC_ surface_elements, DBSPH_find_close_faces, semi_particle_volumes,
!                                         adjacent_faces_isolated_points, area_triangle, DBSPH_kinematics, DBSPH_inlet_outlet, wavy_inlet
!                                      Subroutines/features with relevant modifications:
!                                         SetParticles, area_quadrilateral, Euler, Loop_Irre_2D, Loop_Irre_3D, CalcVarLength, Smoothing p, mixture_viscosity, 
!                                         Gradients_to_MUSCL_boundary, wall_elements_pp
!                                      New input:
!                                         DBSPH (MUSCL_boundary_flag; in_built_monitors; n_monitor_points; n_monitor_regions; n_kinematics_records; dx_dxw;
!                                         k_w; monitor_region(6); monitor_IDs; kinematics(n_records,4); n_inlet; n_outlet; inlet_sections(n_inlet,10): 
!                                         pos(3),normal(3),vel(3),length; outlet_sections(n_outlet,8): pos(3),normal(3), length, pres; 
!                                         surface mesh: surface_mesh_list.txt and surface mesh files; fictitious air reservoirs; parallelepiped domain;
!                                         TyParticle new elements: DBSPH_inlet_ID, DBSPH_outlet_ID
!                                      New output:
!                                         New PV output for DBSPH surface elements: shape coefficient (k_d) and volumes.
!                                         New .txt output file for the global Fx in the selected region (…wall_Fx…txt) 
!                                         New .txt output file for the surface element values in the selected region (…wall_region…txt) 
!                                         New .txt output file for the surface element values, provided selected IDs (…wall_IDs…txt) 
!AA601 end
!
!************************************************************************************
! Module purpose : Main program for calculation flow management
!
! Calling routine: none
!
! Called routines: diagnostic
!                  check_files
!                  Gest_Dealloc
!                  Gest_Input
!                  Gest_Trans
!                  getarg
!                  s_ctime
!                  start_and_stop
!
!************************************************************************************
!
program sphera
!
!.. global modules
!
use FILES_ENTITIES
use GLOBAL_MODULE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE
!
! use OMP_LIB
!
!.. declarations
!
implicit none

character(len=255)     :: nomearg
integer(4)             :: ier, i,n   !standardf90 ,leng,status
integer(4)             :: narg
double precision       :: starttime, endtime
integer(4),external    :: iargc
character(80),external :: lcase
character(len=lencard) :: nomsub = "Sphera"
double precision,external :: omp_get_wtime
!
!.. OMP initializations
!
! integer(4)         :: numthreads
!  numthreads = 2
!  call OMP_SET_NUM_THREADS (numthreads)
!  numthreads = 0
!  numthreads = OMP_GET_NUM_THREADS()
!  write (*,*) 'There are ', numthreads, ' threads'
!stop 
!
!.. Initializations
!
!
!inizializzazione per compatibilita xlf90   icoordp(0:3,2) = (/ 0, 1,3,0,  0, 1,2,3 /)
 icoordp(0,1) = 0
 icoordp(1,1) = 1
 icoordp(2,1) = 3
 icoordp(3,1) = 0
 icoordp(0,2) = 0
 icoordp(1,2) = 1
 icoordp(2,2) = 2
 icoordp(3,2) = 3
!fine inizializzazione per compatibilita xlf90

 pinttimeratio = 0
!
!AA405 rm
! NumPartperLine= 0
! NumPartFace   = 0
!
 NumOpenSides  = 0
 SourceSide    = 0
 Openside      = 0
 NumOpenFaces  = 0
 SourceFace    = 0
 OpenFace      = 0
 mat           = 0
!
!AA406
 itime_jet = 0
! 
! Ncbs          = 0
!
!AA405 sub
 do i=1,MAXOPENSIDES
    RowVelocity(i) = zero
    NumPartperLine(i) = 0
    NumPartFace(i) = 0
 end do
!
 RowPeriod     = zero
 yfila         = zero
 zfila         = zero
 ParticleVolume= zero
 kill_flag = .False.
!
!.. initialization of PgZero
!
 PgZero%vel_type     = "   "
 PgZero%slip         = " "
 PgZero%state        = "   "
 PgZero%kodvel       = 0
 PgZero%koddens      = 0
!!! PgZero%npar2h       = 0
 PgZero%CloseBcOut   = 0
 PgZero%cella        = 0
 PgZero%izona        = 0
 PgZero%icol         = 0
 PgZero%imed         = 0
 PgZero%coord(1)     = zero
 PgZero%coord(2)     = zero
 PgZero%coord(3)     = zero
 PgZero%CoordOld(1)  = zero
 PgZero%CoordOld(2)  = zero
 PgZero%CoordOld(3)  = zero
 PgZero%vel(1)       = zero
 PgZero%vel(2)       = zero
 PgZero%vel(3)       = zero
 PgZero%velass(1)    = zero
 PgZero%velass(2)    = zero
 PgZero%velass(3)    = zero
 PgZero%densass      = zero
 PgZero%acc(1)       = zero
 PgZero%acc(2)       = zero
 PgZero%acc(3)       = zero
 PgZero%mass         = zero
 PgZero%dens         = zero
 PgZero%pres         = zero
 PgZero%dden         = zero
 PgZero%uni          = zero
 PgZero%zer(1)       = zero
 PgZero%zer(2)       = zero
 PgZero%zer(3)       = zero
 PgZero%var(1)       = zero
 PgZero%var(2)       = zero
 PgZero%var(3)       = zero
 PgZero%vpres        = zero
 PgZero%vden         = zero
 PgZero%secinv       = zero
 PgZero%dudx         = zero
 PgZero%dudy         = zero
 PgZero%dvdx         = zero
 PgZero%dvdy         = zero
 PgZero%visc         = zero
 PgZero%mu           = zero
 PgZero%vstart(1)    = zero
 PgZero%vstart(2)    = zero
 PgZero%vstart(3)    = zero
 PgZero%velmorr(1)   = zero
 PgZero%velmorr(2)   = zero
 PgZero%velmorr(3)   = zero
 PgZero%tstop        = zero
 PgZero%mno          = zero
 PgZero%ang          = zero
!.. modello diffusione
 PgZero%VolFra       = zero
 PgZero%rhoc         = zero
 PgZero%rhow         = zero
 PgZero%tiroc        = zero
 PgZero%cden         = zero
 PgZero%wden         = zero
 PgZero%diffu        = zero
 PgZero%coefdif      = zero
 PgZero%veldif(1)    = zero
 PgZero%veldif(2)    = zero
 PgZero%veldif(3)    = zero
!.. modello diffusione
!.. Specific Internal Energy 
 PgZero%IntEn        = zero
 PgZero%Envar        = zero
 PgZero%dEdT         = zero
 PgZero%Csound       = zero
!.. Specific Internal Energy
!
!AA402 start
!.. initialization of ts_pgZero
 ts_pgZero%ts_coord(1)     = zero
 ts_pgZero%ts_coord(2)     = zero
 ts_pgZero%ts_coord(3)     = zero
 ts_pgZero%ts_vel(1)       = zero
 ts_pgZero%ts_vel(2)       = zero
 ts_pgZero%ts_vel(3)       = zero
 ts_pgZero%ts_acc(1)       = zero
 ts_pgZero%ts_acc(2)       = zero
 ts_pgZero%ts_acc(3)       = zero
 ts_pgZero%ts_dden         = zero
 ts_pgZero%ts_dens         = zero
 ts_pgZero%ts_var(1)       = zero
 ts_pgZero%ts_var(2)       = zero
 ts_pgZero%ts_var(3)       = zero
 ts_pgZero%ts_IntEn        = zero
 ts_pgZero%ts_dEdT         = zero
!AA402 end
!
!.. max number of particles
 PARTICLEBUFFER = INIPARTICLEBUFFER
!
!.. inizializzazione vettore indice celle per sub CalcVarLength
!
 indicecelle(1,1) =  0 ; indicecelle(1,2) =  0 ; indicecelle(1,3) =  0
 indicecelle(2,1) =  1 ; indicecelle(2,2) =  0 ; indicecelle(2,3) =  0
 indicecelle(3,1) =  1 ; indicecelle(3,2) =  0 ; indicecelle(3,3) =  1
 indicecelle(4,1) =  0 ; indicecelle(4,2) =  0 ; indicecelle(4,3) =  1
 indicecelle(5,1) = -1 ; indicecelle(5,2) =  0 ; indicecelle(5,3) =  1
 indicecelle(6,1) =  1 ; indicecelle(6,2) =  1 ; indicecelle(6,3) =  0
 indicecelle(7,1) =  0 ; indicecelle(7,2) =  1 ; indicecelle(7,3) =  0
 indicecelle(8,1) = -1 ; indicecelle(8,2) =  1 ; indicecelle(8,3) =  0
 indicecelle(9,1) =  1 ; indicecelle(9,2) =  1 ; indicecelle(9,3) =  1
 indicecelle(10,1) =  0 ; indicecelle(10,2) =  1 ; indicecelle(10,3) =  1
 indicecelle(11,1) = -1 ; indicecelle(11,2) =  1 ; indicecelle(11,3) =  1
 indicecelle(12,1) =  1 ; indicecelle(12,2) =  1 ; indicecelle(12,3) = -1
 indicecelle(13,1) =  0 ; indicecelle(13,2) =  1 ; indicecelle(13,3) = -1
 indicecelle(14,1) = -1 ; indicecelle(14,2) =  1 ; indicecelle(14,3) = -1
!
!.. executable statements
!
 nomecaso = "??????"
 narg     = iargc()
!standardf90 narg = command_argument_count ()
!
!.. check for command line arguments
!
 if (narg < 1) call diagnostic (arg1=1)

 do n = 1,narg
!standardf90   call get_command_argument (n, nomearg, leng, status)
!standardf90   if (status .ne. 0) then
!standardf90     write (*,*) 'get_command_argument failed: status = ', status, ' arg = ', n
!standardf90     stop
!standardf90   end if
!
   do i = 1, len(nomearg)
     nomearg(i:i) = " "
   end do
! modifica per compatibilita con xlf90   call getarg ( int2(n),nomearg )
   call getarg ( n,nomearg )
! fine modifica
!
!.. standard defaults are used, so only the case name has been assigned
!
   if ( nomecaso == "??????" ) then
!standardf90     nomecaso = nomearg(1:leng)
     nomecaso = nomearg
     nomefile(0)  = trim(nomecaso)//".inp"
     nomefile(1)  = trim(nomecaso)//".out"
     nomefile(2)  = trim(nomecaso)//".ris"
     nomefile(3)  = trim(nomecaso)//".rst"
     nomefile(4)  = trim(nomecaso)//".plb"
     nomefile(5)  = trim(nomecaso)//".fro"
     nomefilekill = trim(nomecaso)//".kll"
     nomefileerr  = trim(nomecaso)//".err"
   else
!
!.. erroneous case name assignment. Run is stopped
!
     call diagnostic (arg1=1)
   end if
 end do
!
! .. I/O assignments
!
 call check_files
!
 write (nscr,1000) version, date_vers
 write (nout,1000) version, date_vers
1000 format (10x,' '/  &
              10x,'***************************************************'/  &
              10x,'*                                                 *'/  &
              10x,'*                    SPHERA                       *'/  &
              10x,'*         Smoothed Particle Hydrodynamics         *'/  &
              10x,'*     ----------------------------------------    *'/  &
              10x,'*                                                 *'/  &
              10x,'*         release ',a8," <> ",a14,'      *'/  &
              10x,'*                                                 *'/  &
              10x,'***************************************************'/  &
              10x,'(c)2013 Research on Energy Systems - RSE S.p.A.',/  &
              10x,'Environment and Sustainable Development Dept.',//)
!
 call start_and_stop(0,0)
! 
 write ( nout,* )
 write ( nout,* )  " >> The following files have been assigned and checked:"
 write ( nout,* )  "    Input  file      : ",trim(nomefile(0))
 write ( nout,* )  "    Output file      : ",trim(nomefile(1))
 write ( nout,* )  "    Result file      : ",trim(nomefile(2))
 write ( nout,* )  "    Restart file     : ",trim(nomefile(3))
 write ( nout,* )  "    Free surface file: ",trim(nomefile(4))
 write ( nout,* )  "    Front file       : ",trim(nomefile(5))
 write ( nout,* )
!
 
! .. input deck analysis
!
 call start_and_stop(2,2)
 call Gest_Input
 call start_and_stop(3,2)
!
!.. check limits on particle number
!
 if (PARTICLE_NUMBER_LIMIT /= 0) then
   if (nag > PARTICLE_NUMBER_LIMIT) then
     write ( nscr,* )  "  "
     write ( nscr,* )  " >> The number of particle exceed the maximum allowed."
     write ( nscr,* )  "    Increase the dd parameter in the input file"
     write ( nout,* )  "  "
     write ( nout,* )  " >> The number of particle exceed the maximum allowed."
     write ( nout,* )  "    Increase the dd parameter in the input file"
     stop
   end if
 end if
!
! .. transient loop 
!
! calcolo iniziale tempo Elapsed per versione parallela
 starttime = zero
!$ starttime = omp_get_wtime()
!
 call start_and_stop(2,3)
 call Gest_Trans
 call start_and_stop(3,3)
!
!.. deallocation arrays
!
 call Gest_Dealloc ( nomsub )
!
! calcolo finale tempo Elapsed per versione parallela
 endtime = zero
!$ endtime = omp_get_wtime()
 if (endtime /= zero) then
   write ( nout,* )
   write ( nout,* )
   write ( nout,* ) "Parallel work took ",endtime-starttime," seconds"
   write ( nout,* )
   write ( nout,* )
 end if
!
 write ( nout,* )
!
 call start_and_stop(1,0)
!
 write ( nout,* )
 write ( nout,* )
!
 call s_ctime( nout )
!
 write ( nout,* )
 write ( nout,* )
!
! check if file kill exists. if it exists the execution is stopped saving the results.
 if (kill_flag) then
  open(unit=unitkill,file=nomefilekill,form='formatted',status='old')
  !!!!close (unit=unitkill,disp='delete')
  close (unit=unitkill,status='delete')
  write ( nout,* ) acode,version,date_vers,"- Execution Terminated : kill file found by user request."
 else
  write ( nout,* ) acode,version,date_vers,"- Execution Terminated."
 end if
!
 write ( nout,* )
!
! close the files
!
 if ( nout > 0 ) close ( nout )
 if ( nres > 0 ) close ( nres )
! if ( nsav > 0 ) close ( nsav )
 if ( nplb > 0 ) close ( nplb )
 if ( nfro > 0 ) close ( nfro )
 if ( ncpt > 0 ) close ( ncpt )
!
 stop
!
 contains

!===========================================================================================

subroutine check_files
!
!.. checks for existence and correctness of the case files
!
!.. module assignments
!
!.. implicit declarations
implicit none
!
!.. scalar assignments
!integer(4)   :: n, ier
logical fileexist
!
!.. execution statements
!
!.. check the case file existences
!
! do n = 0, maxfiles
!   inquire ( file=nomefile(n), exist=existfile(n) )
! end do
!
!.. assign the output file
!
 open ( nout, file=nomefile(1), status="unknown", form="formatted", iostat = ier )
 if ( ier /= 0 ) call diagnostic (arg1=2,arg2=1)
!
 fileexist = .false.
!.. check for the input file
 inquire ( file=nomefile(0), exist=fileexist )
! if ( .not. existfile(0) ) then
 if ( .not. fileexist ) then
   call diagnostic (arg1=2,arg2=2)
 else
   open ( ninp, file=nomefile(0), status="old", form="formatted", iostat = ier )
   if ( ier /= 0 ) call diagnostic (arg1=2,arg2=3)
 end if
!
return
!
end subroutine check_files
!
end program sphera
!---split




