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
character(100) :: ainp
integer(4),dimension(20)    :: NumberEntities
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
Npoints = NumberEntities(4)
NumberEntities = 0
! Check the compatibility of the input file with the program version 
current_version = .TRUE.
!------------------------
! Statements
!------------------------
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INPUT VERSION",ninp,nout)) return
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
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"INPUT FILE SECTIONS",ninp,nout))  &
      return
   if (ncord>0) write(nout,"(//,1x,a,/)") lcase(ainp) 
   select case(trim(lcase(trim(ainp))))
      case("##### title #####")
         call ReadInputTitle(ainp,comment,nrighe,ier,ninp,nout)
      case("##### restart #####")
         call ReadInputRestart(ainp,comment,nrighe,ier,ninp,nout)
      case("##### domain #####")
         call ReadInputDomain(NumberEntities,ainp,comment,nrighe,ier,ninp,     &
            nout,nscr)
      case("##### vertices #####")
         call ReadInputVertices(NumberEntities,Vertice, ainp,comment,nrighe,   &
            ier,.TRUE.,ninp,nout)
      case("##### lines #####")
         call ReadInputLines(NumberEntities,BoundaryVertex,Tratto,ainp,        &
            comment,nrighe,ier,ninp,nout)
      case("##### faces #####")
        call ReadInputFaces(NumberEntities,ainp,comment,nrighe,ier,.TRUE.,     &
           ninp,nout)
      case("##### geometry file #####")
         call ReadInputExternalFile(NumberEntities,ainp,comment,nrighe,ier,    &
            OnlyTriangle,ninp,nout,ninp2)
      case("##### bed load transport #####")
         call ReadBedLoadTransport(ainp,comment,nrighe,ier,ninp,nout,nscr)
      case("##### medium #####")
         call ReadInputMedium(NumberEntities,Med,ainp,comment,nrighe,ier,ninp, &
            nout,nscr)
      case("##### body dynamics #####")
         call ReadBodyDynamics(ainp,comment,nrighe,ier,ninp,nout)
! Lower case letters are required
      case("##### dbsph #####") 
         call ReadDBSPH(ainp,comment,nrighe,ier,ninp,nout)
      case("##### boundaries #####")
         call ReadInputBoundaries(NumberEntities,Partz,Tratto,BoundaryVertex,  &
            ainp,comment,nrighe,ier,ninp,nout)
      case("##### run parameters #####")
         call ReadInputRunParameters (ainp,comment,nrighe,ier,ninp,nout,nscr)
      case("##### general physical properties #####")
         call ReadInputGeneralPhysical(NumberEntities,ainp,comment,nrighe,ier, &
            ninp,nout)
      case("##### output regulation #####")
         call ReadInputOutputRegulation (Med,ainp,comment,nrighe,ier,ninp,nout)
      case("##### control points #####")
         call ReadInputControlPoints(NumberEntities,Control_Points,ainp,       &
            comment,nrighe,ier,ninp,nout)
      case("##### control lines #####")
         call ReadInputControlLines(NumberEntities,Control_Points,             &
            Control_Lines,ainp,comment,nrighe,ier, ninp,nout)
      case("##### control sections #####")
         call ReadInputControlSections(NumberEntities,Control_Sections,ainp,   &
            comment,nrighe,ier,ninp,nout)
      case("##### section flow rate #####")
         call ReadSectionFlowRate(ainp,comment,nrighe,ier,ninp,nout)
      case("##### draw options #####")
         call ReadInputDrawOptions(ainp,comment,nrighe,ier,ninp,nout)
      case default 
        ier = 1
   endselect
! Reading error was detected
   if (ier/=0) then
      write(nscr,"(/,1x,a,i8,//)")                                             &
         ">> Reading error was detected in INPUT FILE.  error code= ",ier
      write(nout,"(/,1x,a,i8,//)")                                             &
         ">> Reading error was detected in INPUT FILE.  error code= ",ier
      return
      elseif (ncord>0) then
         write(nout,"(1x,a,/)") lcase(ainp)
   endif
enddo SECTION_LOOP 
! To assign the kernel support 
if (ncord>0) then
   write(nout,"(/,1x,a,//)") ">> END OF INPUT FILE"
   Domain%h = Domain%dx * Domain%trunc
   doubleh = two * Domain%h
   squareh = Domain%h * Domain%h
   doublesquareh = doubleh * doubleh
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

