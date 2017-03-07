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
! Program unit: ReadInputControlSections                      
! Description: Reading control sections (not valid for the flow rate)                       
!-------------------------------------------------------------------------------
subroutine ReadInputControlSections(NumberEntities,Control_Sections,ainp,      &
                                    comment,nrighe,ier,ninp,nout)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,ninp,nout,npts
integer(4),dimension(20) :: NumberEntities
type (TySection),dimension(0:Nsections+1) :: Control_Sections
character(1) :: comment
character(100) :: ainp
integer(4) :: icord,icor2,icor3,icolor,ndiv,ioerr
double precision,dimension(3) :: vp
character(8) :: label
character(100) :: token
character(1),dimension(3) :: CoordLabel = (/ "x", "y", "z" /)
double precision,dimension(3,2) :: values
logical,external :: ReadCheck
integer(4),external :: NumberSectionPoints
character(100),external :: lcase, GetToken
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
if (.NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL SECTIONS DATA",ninp,nout))  &
   return
npts = npoints+npointsl
do while (TRIM(lcase(ainp))/="##### end control sections #####")
   NumberEntities(12) = NumberEntities(12) + 1
   label(1:8) = ainp(1:8)
   vp = zero
   values(:,1) = -99999999.
   values(:,2) =  99999999.
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   token = lcase(GetToken(ainp,1,ioerr))
   if (ioerr==0) then
      select case (token(1:1))
         case ("x","y","z")
            icord = index("xyz",token(1:1))
            token = lcase(GetToken(ainp,2,ioerr))
            if (ioerr==0) read(token,*,iostat=ioerr) vp(icord)
         case default
! Error for unrecognized case 
            ioerr = -1 
      endselect
   endif
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL SECTION "//label//" - CONSTANT COORD. DEFINITION",ninp,nout))   &
      return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   token = lcase(GetToken(ainp,1,ioerr))
   if (ioerr==0) then
      select case (token(1:1))
         case ("x","y","z")
            icor2 = index("xyz",token(1:1))
            if (icor2/=icord) then
               token = lcase(GetToken(ainp,2,ioerr))
               if (ioerr==0) read(token,*,iostat=ioerr) values(icor2,1)
               if (ioerr==0) token = lcase(GetToken(ainp,3,ioerr))
               if (ioerr==0) read(token,*,iostat=ioerr) values(icor2,2)
               else
! Error for unrecognized case 
                  ioerr = -1 
            endif
         case default
! Error for unrecognized case 
            ioerr = -1 
      endselect
   endif
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL SECTION "//label//" - FIRST LIMIT DEFINITION",ninp,nout)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   token = lcase(GetToken(ainp,1,ioerr))
   if (ioerr==0) then
      select case (token(1:1))
         case ("x","y","z")
            icor3 = index("xyz",token(1:1))
            if (icor3/=icord .AND. icor3/=icor2) then
               token = lcase(GetToken(ainp,2,ioerr))
               if (ioerr==0) read(token,*,iostat=ioerr) values(icor3,1)
               if (ioerr==0) token = lcase(GetToken(ainp,3,ioerr))
               if (ioerr==0) read(token,*,iostat=ioerr) values(icor3,2)
               else
! Error for unrecognized case 
                  ioerr = -1 
            endif
         case default
! Error for unrecognized case 
            ioerr = -1 
      endselect
   endif
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL SECTION "//label//" - SECOND LIMIT DEFINITION",ninp,nout))      &
      return
! Colour code: "ptcl"
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   token = lcase(GetToken(ainp,1,ioerr))
   if (ioerr==0) then
      select case (token(1:5))
         case ("color")
            token = lcase(GetToken(ainp,2,ioerr))
            if (ioerr==0) read(token,*,iostat=ioerr) icolor
         case default
! Error for unrecognized case 
            ioerr = -1 
      endselect
   endif
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                                   &
      "CONTROL SECTION "//label//" - COLOR INDEX",ninp,nout)) return
! To estimate the number of points to be generated 
   Ndiv = NumberSectionPoints (values, CoordLabel(icord))
   NumberEntities(13) = NumberEntities(13) + Ndiv
   if (ncord>0) then
      Control_Sections(NumberEntities(12))%Label = label
      Control_Sections(NumberEntities(12))%Tipo = CoordLabel(icord)  
      Control_Sections(NumberEntities(12))%Constant(:) = vp(:)
      Control_Sections(NumberEntities(12))%icont(1) = npts + 1
      Control_Sections(NumberEntities(12))%icont(2) = npts + Ndiv
      Control_Sections(NumberEntities(12))%XYZRange(:,1) = -99999999.
      Control_Sections(NumberEntities(12))%XYZRange(:,2) =  99999999.
      Control_Sections(NumberEntities(12))%XYZRange(icor2,1) =                 &
         minval(values(icor2,:))
      Control_Sections(NumberEntities(12))%XYZRange(icor2,2) =                 &
         maxval(values(icor2,:))
      Control_Sections(NumberEntities(12))%XYZRange(icor3,1) =                 &
         minval(values(icor3,:))
      Control_Sections(NumberEntities(12))%XYZRange(icor3,2) =                 &
         maxval(values(icor3,:)) 
      Control_Sections(NumberEntities(12))%TGLsection = zero
      Control_Sections(NumberEntities(12))%TGLsection(1,1) = one
      Control_Sections(NumberEntities(12))%TGLsection(2,2) = one
      Control_Sections(NumberEntities(12))%TGLsection(3,3) = one
      Control_Sections(NumberEntities(12))%ColorCode = icolor
      call CreateSectionPoints(vp,values,CoordLabel(icord),NumberEntities(12))
      npts = npts + Ndiv
      if (nout>0) then
         write (nout,"(1x,a,i3,1x,a)") "Control section   ",                   &
            NumberEntities(12),                                                &
            "("//Control_Sections(NumberEntities(12))%label//")"
         write (nout,"(1x,a,1x,a,f12.4)") "Constant Coordinate   ",            &
            Control_Sections(NumberEntities(12))%Tipo//"=",                    &
            Control_Sections(NumberEntities(12))%Constant(icord)
         write (nout,"(1x,a,2f12.4)")                                          &
            CoordLabel(icor2)//" Coordinate Limits  ",                         &
            Control_Sections(NumberEntities(12))%XYZRange(icor2,:)
         write (nout,"(1x,a,2f12.4)")                                          &
            CoordLabel(icor3)//" Coordinate Limits  ",                         &
            Control_Sections(NumberEntities(12))%XYZRange(icor3,:) 
         write (nout,"(1x,a,i12)") "First Point:      ",                       &
            Control_Sections(NumberEntities(12))%icont(1)
         write (nout,"(1x,a,i12)") "Last  Point:      ",                       &
            Control_Sections(NumberEntities(12))%icont(2)
      endif
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck (ioerr,ier,nrighe,ainp,"CONTROL SECTIONS DATA",ninp,     &
      nout)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputControlSections

