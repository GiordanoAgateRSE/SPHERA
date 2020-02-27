!-------------------------------------------------------------------------------
! SPHERA v.9.0.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2019 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
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
! Program unit: ReadInput                    
! Description: Reading input data.                    
!-------------------------------------------------------------------------------
subroutine ReadInput(NumberEntities,OnlyTriangle,ier,ainp)
!------------------------
! Modules
!------------------------
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: OnlyTriangle
integer(4) :: ier
character(LEN=lencard) :: ainp
integer(4),dimension(20) :: NumberEntities
logical :: restartOK
integer(4) :: ioerr,nrighe,ioutpo2,iplot_fr,imemo_fr,irest_fr,icpoi_fr,ipllb_fr
integer(4) :: ipllb_md,ioutopt
double precision :: plot_fr,memo_fr,rest_fr,cpoi_fr,pllb_fr
character(1) :: comment = "!"
character(100),external :: lcase,GetToken
logical,external :: ReadCheck
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
restartOK = .FALSE.
ier = 0
vtkconv = .FALSE.
pesodt = zero
dt_alfa_Mon = .false.
ioerr = 0
nrighe = 0
npoints = NumberEntities(4)
NumberEntities = 0
! Check the compatibility of the input file with the program version 
current_version = .TRUE.
!------------------------
! Statements
!------------------------
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INPUT VERSION",ninp,ulog)) return
if (trim(ainp)/=trim(version)) then
   ier = 2
   current_version = .FALSE.
   return
endif
! Loop over "input sections" 
SECTION_LOOP: do while (ioerr==0)
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
! If EOF is reached, then exit, otherwise to check the error code
   if (ioerr==-1) cycle SECTION_LOOP
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INPUT FILE SECTIONS",ninp,ulog))  &
      return
   if (ncord>0) write(ulog,"(//,1x,a,/)") lcase(ainp) 
   select case (trim(lcase(trim(ainp))))
      case("##### title #####")
         call ReadInputTitle(ainp,comment,nrighe,ier,ninp,ulog)
      case("##### restart #####")
         call ReadInputRestart(ainp,comment,nrighe,ier,ninp,ulog)
      case("##### domain #####")
         call ReadInputDomain(NumberEntities,ainp,comment,nrighe,ier,ninp,ulog)
      case("##### vertices #####")
         call ReadInputVertices(NumberEntities,Vertice, ainp,comment,nrighe,   &
            ier,.TRUE.,ninp,ulog)
      case("##### lines #####")
         call ReadInputLines(NumberEntities,BoundaryVertex,Tratto,ainp,        &
            comment,nrighe,ier,ninp,ulog)
      case("##### faces #####")
        call ReadInputFaces(NumberEntities,ainp,comment,nrighe,ier,.TRUE.,     &
           ninp,ulog)
      case("##### geometry file #####")
         call ReadInputExternalFile(NumberEntities,ainp,comment,nrighe,ier,    &
            OnlyTriangle,ninp,ulog,ninp2)
      case("##### bed load transport #####")
         call ReadBedLoadTransport(ainp,comment,nrighe,ier,ninp,ulog,uerr)
      case("##### medium #####")
         call ReadInputMedium(NumberEntities,Med,ainp,comment,nrighe,ier,ninp, &
            ulog)
      case("##### body dynamics #####")
         call ReadBodyDynamics(ainp,comment,nrighe,ier,ninp,ulog)
! Lower case letters are required
      case("##### dbsph #####") 
         call ReadDBSPH(ainp,comment,nrighe,ier,ninp,ulog)
      case("##### boundaries #####")
         call ReadInputBoundaries(NumberEntities,Partz,Tratto,BoundaryVertex,  &
            ainp,comment,nrighe,ier,ninp,ulog)
      case("##### run parameters #####")
         call ReadInputRunParameters (ainp,comment,nrighe,ier,ninp,ulog,uerr)
      case("##### general physical properties #####")
         call ReadInputGeneralPhysical(NumberEntities,ainp,comment,nrighe,ier, &
            ninp,ulog)
      case("##### output regulation #####")
         call ReadInputOutputRegulation (Med,ainp,comment,nrighe,ier,ninp,ulog)
      case("##### control points #####")
         call ReadInputControlPoints(NumberEntities,Control_Points,ainp,       &
            comment,nrighe,ier,ninp,ulog)
      case("##### control lines #####")
         call ReadInputControlLines(NumberEntities,Control_Points,             &
            Control_Lines,ainp,comment,nrighe,ier, ninp,ulog)
      case("##### section flow rate #####")
         call ReadSectionFlowRate(ainp,comment,nrighe,ier)
      case("##### substations #####")
         call ReadSubstations(ainp,comment,nrighe,ier)
      case("##### draw options #####")
         call ReadInputDrawOptions(ainp,comment,nrighe,ier,ninp,ulog)
      case default 
        ier = 1
   endselect
! Reading error was detected
   if (ier/=0) then
      write(uerr,"(/,1x,a,i8,//)")                                             &
         ">> Reading error was detected in INPUT FILE.  error code= ",ier
      write(ulog,"(/,1x,a,i8,//)")                                             &
         ">> Reading error was detected in INPUT FILE.  error code= ",ier
      return
      elseif (ncord>0) then
         write(ulog,"(1x,a,/)") lcase(ainp)
   endif
enddo SECTION_LOOP 
! To assign the kernel support 
if (ncord>0) then
   write(ulog,"(/,1x,a,//)") ">> END OF INPUT FILE"
   Domain%h = Domain%dx * Domain%trunc
   doubleh = two * Domain%h
   squareh = Domain%h * Domain%h
   square_doubleh = doubleh * doubleh
   cubich = squareh * Domain%h
   unosuh = one / Domain%h
   unosusquareh = one / squareh
   eta = 0.001d0 * Domain%h
   eta2 = 0.01d0 * squareh
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInput

