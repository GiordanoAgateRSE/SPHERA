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
! Program unit: ReadInputDrawOptions                       
! Description:                        
!-------------------------------------------------------------------------------
subroutine ReadInputDrawOptions(ainp,comment,nrighe,ier,ninp,ulog)
!------------------------
! Modules
!------------------------
use Static_allocation_module                            
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier, ninp,ulog
character(1) :: comment
character(LEN=lencard) :: ainp
character(4) :: steptime
integer(4) :: ioerr
character(100) :: token
character(100),external :: lcase, GetToken
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
!------------------------
! Statements
!------------------------
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DRAW OPTIONS DATA",ninp,ulog)) return
do while (trim(lcase(ainp))/="##### end draw options #####")
   select case (lcase(GetToken(ainp,1,ioerr)))
      case("vtkconverter")
         token = lcase(GetToken(ainp,(2),ioerr))
         select case (token)
            case("any")
               token = lcase(GetToken(ainp,(3),ioerr))
               read (token,*,iostat=ioerr) freq_time
               if ((ncord>0).and.(ulog>0)) write(ulog,"(1x,a,1pe12.4,a)")      &
                  "VTKConversion any :",freq_time," seconds."
               val_time  = zero  
            case("at")
               token = lcase(GetToken(ainp,(3),ioerr))
               read (token,*,iostat=ioerr) freq_time
               if ((ncord>0).and.(ulog>0)) write(ulog,"(1x,a,1pe12.4,a)")      &
                  "VTKConversion at :",freq_time," second."
                  val_time  = freq_time
                  freq_time = -freq_time
            case("all")
               token = lcase(GetToken(ainp,(3),ioerr))
               read (token,*,iostat=ioerr) steptime
               if (steptime=='time') then
                  freq_time = Domain%memo_fr
                  val_time  = zero  
                  if ((ncord>0).and.(ulog>0)) write(ulog,"(1x,a,1pe12.4,a)")   &
                     "VTKConversion every :",freq_time," second."
                  elseif (steptime=='step') then
                     freq_time = zero
                     val_time  = const_m_9999
                     if ((ncord>0).and.(ulog>0)) write(ulog,"(1x,a)")          &
                        "VTKConversion all steps."
                     else
                        freq_time = Domain%memo_fr
                        val_time  = zero 
                        if ((ncord>0).and.(ulog>0)) write(ulog,                &
                           "(1x,a,1pe12.4,a)") "VTKConversion every :",        &
                           freq_time," second."
               endif
            case default
               freq_time = Domain%memo_fr
               val_time = zero
               if ((ncord)>0.and.(ulog>0)) write(ulog,"(1x,a,1pe12.4,a)")      &
                  "VTKConversion every :",freq_time," second."
         endselect
         if ((ncord>0).and.(ulog>0)) write(ulog,"(1x,a)") " "
         vtkconv = .TRUE.
      case default
         ier = 4
         return
   endselect
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DRAW OPTIONS DATA",ninp,ulog))    &
      return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputDrawOptions

