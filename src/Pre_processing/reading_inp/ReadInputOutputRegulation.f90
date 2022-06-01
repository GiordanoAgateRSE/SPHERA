!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2021 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ReadInputOutputRegulation                          
! Description:                        
!-------------------------------------------------------------------------------

subroutine ReadInputOutputRegulation(Med,ainp,comment,nrighe,ier,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module                                      
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
type (TyMedium), dimension(NMedium) :: Med
integer(4) :: nrighe,ier, ninp,ulog
character(1) :: comment
character(len=lencard) :: ainp
integer(4) :: iplot_fr,imemo_fr,irest_fr,icpoi_fr,ipllb_fr,ipllb_md,n,ioutopt
integer(4) :: ioutpo2,ioerr
double precision :: plot_fr,memo_fr,rest_fr,cpoi_fr,pllb_fr,depth_dt_out
character(100) :: token
character(7),dimension(3) :: outopt = (/ "full   ","partial","null   " /)
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
iplot_fr = 0
plot_fr = zero
imemo_fr = 0
memo_fr = zero
irest_fr = 0
rest_fr = zero
icpoi_fr = 0
cpoi_fr = zero
ipllb_fr = 0
ipllb_md = 0
pllb_fr = zero
ioutopt = 0
ioutpo2 = 3
depth_dt_out = 0.d0
!------------------------
! Statements
!------------------------
call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"OUTPUT REGULATION DATA",ninp,ulog))  &
   return
do while (trim(lcase(ainp))/="##### end output regulation #####")
   select case (trim(lcase(GetToken(1,ainp,ioerr))))
      case ("display")
         token = lcase(GetToken(2,ainp,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "DISPLAY FREQUENCY STEP/TIME",ninp,ulog)) return
         if (token(1:4)=="step") then 
            token = lcase(GetToken(3,ainp,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "DISPLAY FREQUENCY STEP value",ninp,ulog)) return
            read(token,*,iostat=ioerr) iplot_fr
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "DISPLAY FREQUENCY STEP value",ninp,ulog)) return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(3,ainp,ioerr))
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "DISPLAY FREQUENCY TIME value",ninp,ulog)) return
               read(token,*,iostat=ioerr)  plot_fr
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "DISPLAY FREQUENCY TIME value",ninp,ulog)) return
         endif
      case ("results")
         token = lcase(GetToken(2,ainp,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "RESULTS SAVING FREQUENCY STEP/TIME",ninp,ulog)) return
         if (token(1:4)=="step") then 
            token = lcase(GetToken(3,ainp,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESULTS SAVING FREQUENCY STEP value",ninp,ulog)) return
            read(token,*,iostat=ioerr) imemo_fr
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESULTS SAVING FREQUENCY STEP value",ninp,ulog)) return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(3,ainp,ioerr))
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESULTS SAVING FREQUENCY TIME value",ninp,ulog)) return
               read(token,*,iostat=ioerr)  memo_fr
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESULTS SAVING FREQUENCY TIME value",ninp,ulog)) return
         endif
      case ("restart")
         token = lcase(GetToken(2,ainp,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "RESTART SAVING FREQUENCY STEP/TIME",ninp,ulog)) return
         if(token(1:4)=="step") then 
            token = lcase(GetToken(3,ainp,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESTART SAVING FREQUENCY STEP value",ninp,ulog)) return
            read(token,*,iostat=ioerr) irest_fr
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESTART SAVING FREQUENCY STEP value",ninp,ulog)) return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(3,ainp,ioerr))
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESTART SAVING FREQUENCY TIME value",ninp,ulog)) return
               read(token,*,iostat=ioerr)  rest_fr
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESTART SAVING FREQUENCY TIME value",ninp,ulog)) return
         endif
      case ("control")
         token = lcase(GetToken(2,ainp,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "CONTROL POINTS SAVING FREQUENCY STEP/TIME",ninp,ulog)) return
         if(token(1:4)=="step") then 
            token = lcase(GetToken(3,ainp,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "CONTROL POINTS SAVING FREQUENCY STEP value",ninp,ulog))        &
               return
            read(token,*,iostat=ioerr) icpoi_fr
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "CONTROL POINTS SAVING FREQUENCY STEP value",ninp,ulog))        &
               return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(3,ainp,ioerr))
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "CONTROL POINTS SAVING FREQUENCY TIME value",ninp,ulog))     &
                  return
               read(token,*,iostat=ioerr)  cpoi_fr
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "CONTROL POINTS SAVING FREQUENCY TIME value",ninp,ulog))     &
                  return
         endif
      case ("level")
         do n=1,2
            token = lcase(GetToken(2*n,ainp,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL OPTION",ninp,      &
               ulog)) return
            if (token(1:4)=="step") then 
               token = lcase(GetToken(2*n+1,ainp,ioerr))
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "LEVEL STEP FREQUENCY SAVING",ninp,ulog)) return
               read(token,*,iostat=ioerr) ipllb_fr
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "LEVEL STEP FREQUENCY SAVING",ninp,ulog)) return
               elseif (token(1:4)=="time") then 
                  token = lcase(GetToken(2*n+1,ainp,ioerr))
                  if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "LEVEL TIME FREQUENCY SAVING",ninp,ulog)) return
                  read(token,*,iostat=ioerr)  pllb_fr
                  if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "LEVEL TIME FREQUENCY SAVING",ninp,ulog)) return
                  elseif (token(1:6)=="medium") then 
                     token = lcase(GetToken(2*n+1,ainp,ioerr))
                     if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                 &
                        "LEVEL MEDIUM INDEX",ninp,ulog)) return
                     read(token,*,iostat=ioerr) ipllb_md
                     if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                 &
                        "LEVEL MEDIUM INDEX",ninp,ulog)) return
                     if (ipllb_md==0) then
                        ioerr = 5
                        if (.not.ReadCheck(ioerr,ier,nrighe,ainp,              &
                           "LEVEL MEDIUM INDEX",ninp,ulog)) return
                     endif
                     else 
                        ier = 5
                        return
            endif
         enddo
      case ("depth")
         token = lcase(GetToken(2,ainp,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH OPTION",ninp,ulog))   &
            return
         if (token(1:4)=="dt_o") then 
            token = lcase(GetToken(2+1,ainp,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH DT_OUT SAVING",    &
               ninp,ulog)) return
            read(token,*,iostat=ioerr) depth_dt_out
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH DT_OUT SAVING",    &
               ninp,ulog)) return
         endif
      case ("print")
         token = lcase(GetToken(2,ainp,ioerr))
         if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "PRINT MODE AND FREQUENCY FULL/MINIMUM/NULL",ninp,ulog)) return
         if(token(1:4)=="full") then 
            token = lcase(GetToken(3,ainp,ioerr))
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "PRINT MINIMUM FREQUENCY step",ninp,ulog)) return
            read(token,*,iostat=ioerr)  ioutopt
            if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "PRINT MINIMUM FREQUENCY step",ninp,ulog)) return
            ioutopt =-1*abs(ioutopt)
            ioutpo2 = 1
            elseif (token(1:7)=="partial") then 
               token = lcase(GetToken(3,ainp,ioerr))
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "PRINT MINIMUM FREQUENCY step",ninp,ulog)) return
               read(token,*,iostat=ioerr)  ioutopt
               if (.not.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "PRINT MINIMUM FREQUENCY step",ninp,ulog)) return
               ioutpo2 = 2
               elseif (token(1:4)=="null") then 
                  ioutopt = 0
                  ioutpo2 = 3
         endif
      case default
        ier = 4
        return
   endselect
   call ReadRiga(ninp,ainp,ioerr,comment_sym=comment,lines_treated=nrighe)
   if (.not.ReadCheck(ioerr,ier,nrighe,ainp,"OUTPUT REGULATION DATA",ninp,     &
      ulog)) return
enddo
if (input_second_read.eqv..true.) then
   Domain%iplot_fr = iplot_fr
   input_any_t%plot_fr = plot_fr
   Domain%imemo_fr = imemo_fr
   input_any_t%memo_fr = memo_fr
   Domain%irest_fr = irest_fr
   input_any_t%rest_fr = rest_fr
   Domain%icpoi_fr = icpoi_fr
   input_any_t%cpoi_fr = cpoi_fr
   input_any_t%ipllb_md = ipllb_md
   Domain%ipllb_fr = ipllb_fr
   input_any_t%pllb_fr = pllb_fr
   input_any_t%depth_dt_out = depth_dt_out
   Domain%depth_it_out_last = 0    
   Domain%ioutopt  = ioutopt
   if (ulog>0) then 
      write(ulog,"(1x,a,i12)")   "Displaying     step frequency : ",           &
         Domain%iplot_fr
      write(ulog,"(1x,a,e12.4)") "Displaying     time frequency : ",           &
         input_any_t%plot_fr
      write(ulog,"(1x,a,i12)")   "Results saving step frequency : ",           &
         Domain%imemo_fr
      write(ulog,"(1x,a,e12.4)") "Results saving time frequency : ",           &
         input_any_t%memo_fr
      write(ulog,"(1x,a,i12)")   "Restart saving step frequency : ",           &
         Domain%irest_fr
      write(ulog,"(1x,a,e12.4)") "Restart saving time frequency : ",           &
         input_any_t%rest_fr
      write(ulog,"(1x,a,i12)")   "Ctrl.Points    step frequency : ",           &
         Domain%icpoi_fr
      write(ulog,"(1x,a,e12.4)") "Ctrl.Points    time frequency : ",           &
         input_any_t%cpoi_fr
      write(ulog,"(1x,a,i12)")   "Level Medium Index            : ",           &
         input_any_t%ipllb_md
      if (input_any_t%ipllb_md>0) write(ulog,"(1x,a,e12.4)")                   &
                                 "Level Medium Density Limit    : ",           &
                                 med(input_any_t%ipllb_md)%den0 * half
      write(ulog,"(1x,a,i12)")   "Level          step frequency : ",           &
         Domain%ipllb_fr
      write(ulog,"(1x,a,e12.4)") "Level          time frequency : ",           &
         input_any_t%pllb_fr
      write(ulog,"(1x,a,e12.4)") "Depth          dt_out(s)      : ",           &
         input_any_t%depth_dt_out
      write(ulog,"(1x,a,a)")     "Printing Option               : ",           &
         outopt(ioutpo2)
      if (Domain%ioutopt>0) write(ulog,"(1x,a,i12)")                           &
         "Printing       step frequency : ",iabs(Domain%ioutopt)
      write(ulog,"(1x,a)") " "
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputOutputRegulation
