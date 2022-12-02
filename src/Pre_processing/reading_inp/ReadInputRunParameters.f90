!-------------------------------------------------------------------------------
! SPHERA v.10.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2022 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
! This file is part of SPHERA v.10.0.0
! SPHERA v.10.0.0 is free software: you can redistribute it and/or modify
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
! Program unit: ReadInputRunParameters                            
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputRunParameters(ainp,comment,nrighe,ier,ninp,ulog,uerr)
!------------------------
! Modules
!------------------------
use Static_allocation_module                                   
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: ALE3,C1_BE,C1_monitors,neg_p_FS,C1_divu,C1_ALE1eCE,ALE3_vol_corr
integer(4) :: nrighe,ier,ninp,ulog,uerr
character(1) :: comment
character(len=lencard) :: ainp
integer(4) :: itmax
double precision :: tmax,CFL,TetaP,TetaV,COEFNMAXPARTJ,COEFNMAXPARTI,vsc_coeff
double precision :: invB_det_thresh,max_delta_mass
integer(4) :: ioerr,time_split,RKscheme,body_part_reorder
#ifdef SPACE_3D
integer(4) :: MAXCLOSEBOUNDFACES,MAXNUMCONVEXEDGES,GCBFVecDim_loc,nag_aux
#endif
integer(4) :: density_thresholds
character(100) :: token
logical,external :: ReadCheck
character(100),external :: lcase
!------------------------
! Explicit interfaces
!------------------------
interface
   subroutine ReadRiga(ninp,ainp,io_err,comment_sym,lines_treated)
      implicit none
      integer(4),intent(in) :: ninp
      character(*),intent(inout) :: ainp
      integer(4),intent(out) :: io_err
      character(1),intent(in),optional :: comment_sym
      integer(4),intent(inout),optional :: lines_treated
   end subroutine ReadRiga
   character(100) function GetToken(itok,ainp,io_err)
      implicit none
      integer(4),intent(in) :: itok
      character(*),intent(in) :: ainp
      integer(4),intent(out) :: io_err
   end function GetToken
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"RUN PARAMETERS DATA",ninp,ulog))     &
   return
do while (trim(lcase(ainp))/="##### end run parameters #####")
   read(ainp,*,iostat=ioerr) tmax
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"MAX. TRANSIENT TIME & ITERATIONS",&
      ninp,ulog)) return
   if (ioerr==0) then
      token = GetToken(2,ainp,ioerr)
      if (ioerr==0) then
         read(token,*,iostat=ioerr) itmax
         else
            itmax = 1000000000
            ioerr = 0
      endif
   endif
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) CFL,vsc_coeff,time_split,RKscheme,pesodt,         &
      dt_alfa_Mon
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"TIME INTEGRATION",ninp,ulog))     &
      return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) TetaP,TetaV
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"TETAP & TETAV",ninp,ulog)) return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) COEFNMAXPARTI,COEFNMAXPARTJ,body_part_reorder
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"COEFNMAXPARTI and COEFNMAXPARTJ ",&
      ninp,ulog)) return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) ALE3,C1_BE,C1_divu,C1_ALE1eCE,C1_monitors,        &
      invB_det_thresh,neg_p_FS,max_delta_mass,ALE3_vol_corr
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
"ALE3,C1_BE,C1_divu,C1_ALE1eCE,C1_monitors,invB_det_thresh,neg_p_FS,max_delta_mass,ALE3_vol_corr",&
         ninp,ulog)) return
#ifdef SPACE_3D
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
      read(ainp,*,iostat=ioerr) nag_aux,MAXCLOSEBOUNDFACES,MAXNUMCONVEXEDGES,  &
         GCBFVecDim_loc
      if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                                &
         "NAG_AUX, MAXCLOSEBOUNDFACES, MAXNUMCONVEXEDGES, GCBFVECDIM_LOC ",    &
         ninp,ulog)) return
#endif
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   read(ainp,*,iostat=ioerr) density_thresholds
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DENSITY_THRESHOLDS ",ninp,ulog))  &
      return
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"RUN PARAMETERS DATA",ninp,ulog))  &
      return
enddo
! Assigning the values read
if (input_second_read.eqv..true.) then
   input_any_t%tmax = tmax
   input_any_t%itmax = itmax
   input_any_t%CFL = CFL
   input_any_t%vsc_coeff = vsc_coeff
   Domain%time_split = time_split
   if (time_split==1) then
      Domain%RKscheme = 1
      else
      if (RKscheme==0) then
         write(uerr,"(1x,a)") " "
         write(ulog,"(1x,a)") " "
         write(uerr,"(1x,a)")                                                  &
            "Error setting run parameters. RKscheme must be greather than 0."
         write(ulog,"(1x,a)")                                                  &
            "Error setting run parameters. RKscheme must be greather than 0."
         stop
      endif
      Domain%RKscheme  = RKscheme
   endif
   input_any_t%TetaP = TetaP
   input_any_t%TetaV = TetaV
   Domain%COEFNMAXPARTI = COEFNMAXPARTI
   input_any_t%COEFNMAXPARTJ = COEFNMAXPARTJ
   input_any_t%body_part_reorder = body_part_reorder
   input_any_t%ALE3 = ALE3
   input_any_t%C1_BE = C1_BE
   input_any_t%C1_divu = C1_divu
   input_any_t%C1_ALE1eCE = C1_ALE1eCE
   input_any_t%C1_monitors = C1_monitors
   input_any_t%invB_det_thresh = invB_det_thresh
   input_any_t%neg_p_FS = neg_p_FS
   input_any_t%ALE3_vol_corr = ALE3_vol_corr
   input_any_t%max_delta_mass = max_delta_mass
#ifdef SPACE_3D
      input_any_t%MAXCLOSEBOUNDFACES = MAXCLOSEBOUNDFACES
      input_any_t%MAXNUMCONVEXEDGES = MAXNUMCONVEXEDGES
      if (.not.restart) GCBFVecDim = GCBFVecDim_loc
      Domain%nag_aux = nag_aux
#endif
   input_any_t%density_thresholds = density_thresholds
   if (ulog>0) then
      write(ulog,"(1x,a,1p,e12.4)") "TMAX                       : ",           &
         input_any_t%tmax
      write(ulog,"(1x,a,i12)")      "ITMAX                      : ",           &
         input_any_t%itmax
      write(ulog,"(1x,a,1p,e12.4)") "CFL                        : ",           &
         input_any_t%CFL
      write(ulog,"(1x,a,1p,e12.4)") "vsc_coeff                  : ",           &
         input_any_t%vsc_coeff
      write(ulog,"(1x,a,1p,i1)")    "staggering option          : ",           &
         Domain%time_split
      write(ulog,"(1x,a,1p,i1)")    "RKscheme                   : ",           &
         Domain%RKscheme
      write(ulog,"(1x,a,1p,e12.4)") "SMOOTHING PRES             : ",           &
         input_any_t%TetaP
      write(ulog,"(1x,a,1p,e12.4)") "SMOOTHING VEL              : ",           &
         input_any_t%TetaV
      write(ulog,"(1x,a,1p,e12.4)") "COEFNMAXPARTI              : ",           &
         Domain%COEFNMAXPARTI
      write(ulog,"(1x,a,1p,e12.4)") "COEFNMAXPARTJ              : ",           &
         input_any_t%COEFNMAXPARTJ
      write(ulog,"(1x,a,1p,i1)")    "body_part_reorder          : ",           &
         input_any_t%body_part_reorder
      write(ulog,"(1x,a,1p,l8)")    "ALE3               : ",                   &
         input_any_t%ALE3
      write(ulog,"(1x,a,1p,l8)")    "C1_BE              : ",                   &
         input_any_t%C1_BE
      if (input_any_t%C1_BE) then
         write(ulog,"(1x,a,1p,l8)")    "C1_divu            : ",                &
            input_any_t%C1_divu
         write(ulog,"(1x,a,1p,l8)")    "C1_ALE1eCE         : ",                &
            input_any_t%C1_ALE1eCE
      endif
      write(ulog,"(1x,a,1p,l8)")    "C1_monitors               : ",            &
         input_any_t%C1_monitors
      if (input_any_t%C1_BE) then
         write(ulog,"(1x,a,1p,e12.4)") "invB_det_thresh           : ",         &
            input_any_t%invB_det_thresh
      endif
      if (input_any_t%ALE3) then
         write(ulog,"(1x,a,1p,l8)")    "neg_p_FS                  : ",         &
            input_any_t%neg_p_FS
         write(ulog,"(1x,a,1p,e12.4)")    "max_delta_mass            : ",      &
            input_any_t%max_delta_mass
         write(ulog,"(1x,a,1p,l8)")    "ALE3_vol_corr             : ",         &
            input_any_t%ALE3_vol_corr
      endif
#ifdef SPACE_3D
      write(ulog,"(1x,a,1p,i12)")   "NAG_AUX                    : ",           &
         Domain%nag_aux
      write(ulog,"(1x,a,1p,i12)")   "MAXCLOSEBOUNDFACES         : ",           &
         input_any_t%MAXCLOSEBOUNDFACES
      write(ulog,"(1x,a,1p,i12)")   "MAXNUMCONVEXEDGES          : ",           &
         input_any_t%MAXNUMCONVEXEDGES
      if (.not.restart) then
      write(ulog,"(1x,a,1p,i12)")   "GCBFVecDim (input file)    : ",           &
         GCBFVecDim
      endif
#endif
      write(ulog,"(1x,a,1p,i1)")    "density_thresholds         : ",           &
         input_any_t%density_thresholds       
      write(ulog,"(1x,a)") " "
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputRunParameters
