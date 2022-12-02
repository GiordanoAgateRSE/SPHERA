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
! Program unit: ReadRiga                               
! Description: This program unit reads the first non-commented non-empty line 
!              of the I/O unit "ninp" and assign it to the character variable 
!              "ainp". The comment symbol is "comment_sym". The 
!              "Tab" symbols are replaced with blank spaces. All the symbols 
!              after a comment symbol (included) are replaced with blank spaces.
!              The number of lines treated is "lines_treated", which is 
!              equal to one (the only line read) plus the number of lines 
!              skipped. The comment symbol and the number of lines treated are 
!              optional arguments.
!-------------------------------------------------------------------------------
subroutine ReadRiga(ninp,ainp,io_err,comment_sym,lines_treated)
!------------------------
! Modules
!------------------------
implicit none
integer(4),intent(in) :: ninp
character(*),intent(inout) :: ainp
integer(4),intent(out) :: io_err
character(1),intent(in),optional :: comment_sym
integer(4),intent(inout),optional :: lines_treated
integer(4) :: n,l
!------------------------
! Declarations
!------------------------
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
io_err = 0
!------------------------
! Statements
!------------------------
do while (io_err==0)
   read(ninp,"(a)",iostat=io_err) ainp
   if (present(lines_treated)) lines_treated = lines_treated + 1
   if ((io_err==0).and.(trim(ainp)/="")) then
      if (present(comment_sym)) then
! Skip a comment line
         if (ainp(1:1)==comment_sym) cycle
! Line with a comment after the source-code tokens: all the symbols after a 
! comment symbol (included) are replaced with blank spaces
         l = index(ainp,comment_sym)
         if (l>0) then
            do n=l,len(trim(ainp))
               ainp(n:n) = " "
            enddo
         endif
      endif
! Replacing tabs with blank spaces
      do n=1,len(trim(ainp))
         if (iachar(ainp(n:n))==9) ainp(n:n) = " "
      enddo
      exit
   endif
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadRiga
