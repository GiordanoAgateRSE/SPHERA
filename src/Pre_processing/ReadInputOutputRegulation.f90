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
! Program unit: ReadInputOutputRegulation                          
! Description:                        
!-------------------------------------------------------------------------------

subroutine ReadInputOutputRegulation(Med,ainp,comment,nrighe,ier,ninp,nout)
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
integer(4) :: nrighe,ier, ninp,nout
character(1) :: comment
character(100) :: ainp
integer(4) :: iplot_fr,imemo_fr,irest_fr,icpoi_fr,ipllb_fr,ipllb_md,n,ioutopt
integer(4) :: ioutpo2,ioerr
double precision :: plot_fr,memo_fr,rest_fr,cpoi_fr,pllb_fr,depth_dt_out
character(100) :: token
character(7),dimension(3) :: outopt = (/ "full   ","partial","null   " /)
logical,external :: ReadCheck
character(100),external :: GetToken,lcase
!------------------------
! Explicit interfaces
!------------------------
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
depth_dt_out = 0.0d0
!------------------------
! Statements
!------------------------
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"OUTPUT REGULATION DATA",ninp,nout))  &
   return
do while (TRIM(lcase(ainp))/="##### end output regulation #####")
   select case (TRIM(lcase(GetToken(ainp,1,ioerr))))
      case ("display")
         token = lcase(GetToken(ainp,2,ioerr))
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "DISPLAY FREQUENCY STEP/TIME",ninp,nout)) return
         if (token(1:4)=="step") then 
            token = lcase(GetToken(ainp,3,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "DISPLAY FREQUENCY STEP value",ninp,nout)) return
            read(token,*,iostat=ioerr) iplot_fr
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "DISPLAY FREQUENCY STEP value",ninp,nout)) return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(ainp,3,ioerr))
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "DISPLAY FREQUENCY TIME value",ninp,nout)) return
               read(token,*,iostat=ioerr)  plot_fr
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "DISPLAY FREQUENCY TIME value",ninp,nout)) return
         endif
      case ("results")
         token = lcase(GetToken(ainp,2,ioerr))
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "RESULTS SAVING FREQUENCY STEP/TIME",ninp,nout)) return
         if (token(1:4)=="step") then 
            token = lcase(GetToken(ainp,3,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESULTS SAVING FREQUENCY STEP value",ninp,nout)) return
            read(token,*,iostat=ioerr) imemo_fr
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESULTS SAVING FREQUENCY STEP value",ninp,nout)) return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(ainp,3,ioerr))
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESULTS SAVING FREQUENCY TIME value",ninp,nout)) return
               read(token,*,iostat=ioerr)  memo_fr
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESULTS SAVING FREQUENCY TIME value",ninp,nout)) return
         endif
      case ("restart")
         token = lcase(GetToken(ainp,2,ioerr))
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "RESTART SAVING FREQUENCY STEP/TIME",ninp,nout)) return
         if(token(1:4)=="step") then 
            token = lcase(GetToken(ainp,3,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESTART SAVING FREQUENCY STEP value",ninp,nout)) return
            read(token,*,iostat=ioerr) irest_fr
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "RESTART SAVING FREQUENCY STEP value",ninp,nout)) return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(ainp,3,ioerr))
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESTART SAVING FREQUENCY TIME value",ninp,nout)) return
               read(token,*,iostat=ioerr)  rest_fr
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "RESTART SAVING FREQUENCY TIME value",ninp,nout)) return
         endif
      case ("control")
         token = lcase(GetToken(ainp,2,ioerr))
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "CONTROL POINTS SAVING FREQUENCY STEP/TIME",ninp,nout)) return
         if(token(1:4)=="step") then 
            token = lcase(GetToken(ainp,3,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "CONTROL POINTS SAVING FREQUENCY STEP value",ninp,nout))        &
               return
            read(token,*,iostat=ioerr) icpoi_fr
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "CONTROL POINTS SAVING FREQUENCY STEP value",ninp,nout))        &
               return
            elseif (token(1:4)=="time") then 
               token = lcase(GetToken(ainp,3,ioerr))
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "CONTROL POINTS SAVING FREQUENCY TIME value",ninp,nout))     &
                  return
               read(token,*,iostat=ioerr)  cpoi_fr
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "CONTROL POINTS SAVING FREQUENCY TIME value",ninp,nout))     &
                  return
         endif
      case ("level")
         do n=1,2
            token = lcase(GetToken(ainp,2*n,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"LEVEL OPTION",ninp,      &
               nout)) return
            if (token(1:4)=="step") then 
               token = lcase(GetToken(ainp,2*n+1,ioerr))
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "LEVEL STEP FREQUENCY SAVING",ninp,nout)) return
               read(token,*,iostat=ioerr) ipllb_fr
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "LEVEL STEP FREQUENCY SAVING",ninp,nout)) return
               elseif (token(1:4)=="time") then 
                  token = lcase(GetToken(ainp,2*n+1,ioerr))
                  if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "LEVEL TIME FREQUENCY SAVING",ninp,nout)) return
                  read(token,*,iostat=ioerr)  pllb_fr
                  if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                    &
                     "LEVEL TIME FREQUENCY SAVING",ninp,nout)) return
                  elseif (token(1:6)=="medium") then 
                     token = lcase(GetToken(ainp,2*n+1,ioerr))
                     if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                 &
                        "LEVEL MEDIUM INDEX",ninp,nout)) return
                     read(token,*,iostat=ioerr) ipllb_md
                     if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                 &
                        "LEVEL MEDIUM INDEX",ninp,nout)) return
                     if (ipllb_md==0) then
                        ioerr = 5
                        if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,              &
                           "LEVEL MEDIUM INDEX",ninp,nout)) return
                     endif
                     else 
                        ier = 5
                        return
            endif
         enddo
      case ("depth")
         token = lcase(GetToken(ainp,2,ioerr))
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH OPTION",ninp,nout))   &
            return
         if (token(1:4)=="dt_o") then 
            token = lcase(GetToken(ainp,2+1,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH DT_OUT SAVING",    &
               ninp,nout)) return
            read(token,*,iostat=ioerr) depth_dt_out
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DEPTH DT_OUT SAVING",    &
               ninp,nout)) return
         endif
      case ("print")
         token = lcase(GetToken(ainp,2,ioerr))
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "PRINT MODE AND FREQUENCY FULL/MINIMUM/NULL",ninp,nout)) return
         if(token(1:4)=="full") then 
            token = lcase(GetToken(ainp,3,ioerr))
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "PRINT MINIMUM FREQUENCY step",ninp,nout)) return
            read(token,*,iostat=ioerr)  ioutopt
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "PRINT MINIMUM FREQUENCY step",ninp,nout)) return
            ioutopt =-1*abs(ioutopt)
            ioutpo2 = 1
            elseif (token(1:7)=="partial") then 
               token = lcase(GetToken(ainp,3,ioerr))
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "PRINT MINIMUM FREQUENCY step",ninp,nout)) return
               read(token,*,iostat=ioerr)  ioutopt
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "PRINT MINIMUM FREQUENCY step",ninp,nout)) return
               ioutpo2 = 2
               elseif (token(1:4)=="null") then 
                  ioutopt = 0
                  ioutpo2 = 3
         endif
      case default
        ier = 4
        return
   endselect
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"OUTPUT REGULATION DATA",ninp,     &
      nout)) return
enddo
if (ncord>0) then
   Domain%iplot_fr = iplot_fr
   Domain%plot_fr = plot_fr
   Domain%imemo_fr = imemo_fr
   Domain%memo_fr = memo_fr
   Domain%irest_fr = irest_fr
   Domain%rest_fr = rest_fr
   Domain%icpoi_fr = icpoi_fr
   Domain%cpoi_fr = cpoi_fr
   Domain%ipllb_md = ipllb_md
   Domain%ipllb_fr = ipllb_fr
   Domain%pllb_fr = pllb_fr
   Domain%depth_dt_out = depth_dt_out
   Domain%depth_it_out_last = 0.d0    
   Domain%ioutopt  = ioutopt
   if (nout>0) then 
      write(nout,"(1x,a,i12)")   "Displaying     step frequency : ",           &
         Domain%iplot_fr
      write(nout,"(1x,a,e12.4)") "Displaying     time frequency : ",           &
         Domain%plot_fr
      write(nout,"(1x,a,i12)")   "Results saving step frequency : ",           &
         Domain%imemo_fr
      write(nout,"(1x,a,e12.4)") "Results saving time frequency : ",           &
         Domain%memo_fr
      write(nout,"(1x,a,i12)")   "Restart saving step frequency : ",           &
         Domain%irest_fr
      write(nout,"(1x,a,e12.4)") "Restart saving time frequency : ",           &
         Domain%rest_fr
      write(nout,"(1x,a,i12)")   "Ctrl.Points    step frequency : ",           &
         Domain%icpoi_fr
      write(nout,"(1x,a,e12.4)") "Ctrl.Points    time frequency : ",           &
         Domain%cpoi_fr
      write(nout,"(1x,a,i12)")   "Level Medium Index            : ",           &
         Domain%ipllb_md
      if (Domain%ipllb_md>0) write(nout,"(1x,a,e12.4)")                        &
                                 "Level Medium Density Limit    : ",           &
                                 med(Domain%ipllb_md)%den0 * half
      write(nout,"(1x,a,i12)")   "Level          step frequency : ",           &
         Domain%ipllb_fr
      write(nout,"(1x,a,e12.4)") "Level          time frequency : ",           &
         Domain%pllb_fr
      write(nout,"(1x,a,e12.4)") "Depth          dt_out(s)      : ",           &
         Domain%depth_dt_out
      write(nout,"(1x,a,a)")     "Printing Option               : ",           &
         outopt(ioutpo2)
      if (Domain%ioutopt>0) write(nout,"(1x,a,i12)")                           &
         "Printing       step frequency : ",iabs(Domain%ioutopt)
      write(nout,"(1x,a)") " "
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputOutputRegulation

