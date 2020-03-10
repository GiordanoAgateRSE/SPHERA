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
integer(4) :: nrighe,ier,ninp,ulog,uerr
character(1) :: comment
character(LEN=lencard) :: ainp
integer(4) :: itmax
double precision :: tmax,CFL,TetaP,TetaV,COEFNMAXPARTJ,COEFNMAXPARTI,vsc_coeff
integer(4) :: ioerr,time_split,RKscheme,body_part_reorder,MAXCLOSEBOUNDFACES
integer(4) :: MAXNUMCONVEXEDGES,GCBFVecDim_loc,density_thresholds,nag_aux
character(1) :: Psurf
character(100) :: token
logical,external :: ReadCheck
character(100),external :: lcase,GetToken
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"RUN PARAMETERS DATA",ninp,ulog))     &
   return
do while (trim(lcase(ainp))/="##### end run parameters #####")
   read (ainp,*,iostat=ioerr) tmax
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MAX. TRANSIENT TIME & ITERATIONS",&
      ninp,ulog)) return
   if (ioerr==0) then
      token = GetToken(ainp,2,ioerr)
      if (ioerr==0) then
         read (token,*,iostat=ioerr) itmax
         else
            itmax = 1000000000
            ioerr = 0
      endif
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) CFL,vsc_coeff,time_split,RKscheme,pesodt,        &
      dt_alfa_Mon
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TIME INTEGRATION",ninp,ulog))     &
      return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) TetaP,TetaV
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TETAP & TETAV",ninp,ulog)) return
   token = GetToken(ainp,3,ioerr)
   if (ioerr==0) then
      read (token,*,iostat=ioerr) Psurf
      Psurf = trim(lcase(Psurf))
      if ((Psurf/='o').and.(Psurf/='s').and.(Psurf/='a')) then
         write(uerr,"(1x,a)")                                                  &
            "Error setting run parameters. SMOOTHING Pressure Surface not set."
         write(ulog,"(1x,a)")                                                  &
            "Error setting run parameters. SMOOTHING Pressure Surface not set."
         stop
      endif
      else
         if (ulog>0) write(ulog,*) "Unknown option: ",trim(ainp),              &
            " in run parameters."
         write(uerr,*) "Unknown option: ",trim(ainp)," in run parameters."
         stop
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) COEFNMAXPARTI,COEFNMAXPARTJ,body_part_reorder
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"COEFNMAXPARTI and COEFNMAXPARTJ ",&
      ninp,ulog)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) nag_aux,MAXCLOSEBOUNDFACES,MAXNUMCONVEXEDGES,    &
      GCBFVecDim_loc
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "NAG_AUX, MAXCLOSEBOUNDFACES, MAXNUMCONVEXEDGES, GCBFVECDIM_LOC ",ninp,  &
      ulog)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read (ainp,*,iostat=ioerr) density_thresholds
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DENSITY_THRESHOLDS ",ninp,ulog))  &
      return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"RUN PARAMETERS DATA",ninp,ulog))  &
      return
enddo
! Assigning the values read
if (ncord>0) then
   Domain%tmax = tmax
   Domain%itmax = itmax
   Domain%CFL = CFL
   Domain%vsc_coeff = vsc_coeff
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
   Domain%TetaP = TetaP
   Domain%TetaV = TetaV
   Domain%Psurf = Psurf
   Domain%COEFNMAXPARTI = COEFNMAXPARTI
   Domain%COEFNMAXPARTJ = COEFNMAXPARTJ
   Domain%body_part_reorder = body_part_reorder
   Domain%MAXCLOSEBOUNDFACES = MAXCLOSEBOUNDFACES
   Domain%MAXNUMCONVEXEDGES = MAXNUMCONVEXEDGES
   Domain%nag_aux = nag_aux
   if (.not.restart) GCBFVecDim = GCBFVecDim_loc
   Domain%density_thresholds = density_thresholds
   if (ulog>0) then
      write(ulog,"(1x,a,1p,e12.4)") "TMAX                       : ",           &
         Domain%tmax
      write(ulog,"(1x,a,i12)")      "ITMAX                      : ",           &
         Domain%itmax
      write(ulog,"(1x,a,1p,e12.4)") "CFL                        : ",           &
         Domain%CFL
      write(ulog,"(1x,a,1p,e12.4)") "vsc_coeff                  : ",           &
         Domain%vsc_coeff
      write(ulog,"(1x,a,1p,i1)")    "staggering option          : ",           &
         Domain%time_split
      write(ulog,"(1x,a,1p,i1)")    "RKscheme                   : ",           &
         Domain%RKscheme
      write(ulog,"(1x,a,1p,e12.4)") "SMOOTHING PRES             : ",           &
         Domain%TetaP
      write(ulog,"(1x,a,1p,e12.4)") "SMOOTHING VEL              : ",           &
         Domain%TetaV
      if (Domain%Psurf=='o') then
         write(ulog,"(1x,a,a)") "SMOOTHING Pressure Surface : ","Original"
         elseif (Domain%Psurf=='s') then
            write(ulog,"(1x,a,a)")                                             &
"SMOOTHING Pressure Surface : ","Delta Pressure from hydrostatic"
            elseif (Domain%Psurf=='a') then
               write(ulog,"(1x,a,a)")                                          &
"SMOOTHING Pressure Surface : ","Weight calculation with atmospheric pressure"
      endif
      write(ulog,"(1x,a,1p,e12.4)") "COEFNMAXPARTI              : ",           &
         Domain%COEFNMAXPARTI
      write(ulog,"(1x,a,1p,e12.4)") "COEFNMAXPARTJ              : ",           &
         Domain%COEFNMAXPARTJ
      write(ulog,"(1x,a,1p,i1)")    "body_part_reorder          : ",           &
         Domain%body_part_reorder
      write(ulog,"(1x,a,1p,i12)")   "NAG_AUX                    : ",           &
         Domain%nag_aux
      write(ulog,"(1x,a,1p,i12)")   "MAXCLOSEBOUNDFACES         : ",           &
         Domain%MAXCLOSEBOUNDFACES
      write(ulog,"(1x,a,1p,i12)")   "MAXNUMCONVEXEDGES          : ",           &
         Domain%MAXNUMCONVEXEDGES
      if (.not.restart) then
      write(ulog,"(1x,a,1p,i12)")   "GCBFVecDim (input file)    : ",           &
         GCBFVecDim
      endif
      write(ulog,"(1x,a,1p,i1)")    "density_thresholds         : ",           &
         Domain%density_thresholds       
      write(ulog,"(1x,a)") " "
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputRunParameters

