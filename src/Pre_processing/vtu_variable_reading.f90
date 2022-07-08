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
! Program unit: vtu_variable_reading
! Description: Decoding of a ".vtu" variable within a ".vtu" file. The 
!              recognized ".vtu" variables are: 
!                 var_name="types": ".vtu" cell type;
!                 var_name="connectivity": vertex ID belonging to the ".vtu" 
!                                          cell;
!                 var_name="Points": ".vtu" point position;
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Further Copyright acknowledgments
! The first version of this program unit represented a modification of an 
! extract of the software tool vtu2FEHM (RSE SpA). Please refer to the git log 
! file for further details.
!-------------------------------------------------------------------------------
#if (defined SPACE_3D) && (defined SOLID_BODIES)
subroutine vtu_variable_reading(i_vtu_grid,n_nodes,n_components,n_lines_var,   &
   n_tokens_line,var_name)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_vtu_grid,n_nodes,n_components,n_lines_var
integer(4),intent(in) :: n_tokens_line
character(100),intent(in) :: var_name
integer(4) :: n_values,i_line,io_err,i_tok,i_value,i_node,i_component
character(100) :: token
character(200) :: line
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
n_values = n_components * n_nodes
!------------------------
! Statements
!------------------------
! Loop over all the variable records/lines
do i_line=1,n_lines_var
! Read the line
   call ReadRiga(uvtu,line,io_err)
! Read the tokens, which are here the values of the variable components
   do i_tok=1,n_tokens_line
      token = GetToken(i_tok,line,io_err)
! This assessment of "i_component" and "i_value" seems faster (one "if" 
! contruct less) than equivalent incremental countings
      i_value = (i_line - 1) * n_tokens_line + i_tok
      i_node = int((i_value - 1) / n_components) + 1
      i_component = mod(i_value,n_components)
      if (i_component==0) i_component = n_components
! Check for errors or end of dataset
      if (io_err/=0) then
         if(i_value<=n_values) then
! Reading error 
            write(uerr,*) 'Error in reading the value n. ',i_tok,              &
               'of the line n. ',i_line,'for the current variable ',var_name,  &
               'in the program unit "vtu_variable_reading". The program stops',&
               ' here.'
            stop
            else
! Reading completed. Last record/line: the number of values can be smaller than 
! the previous lines.
               exit
         endif
      endif
      select case (trim(adjustl(var_name)))
         case ("Points")
            read(token,*)                                                      &
               vtu_grids(i_vtu_grid)%points%vertex(i_node)%pos(i_component)
         case ("types")
            read(token,*) vtu_grids(i_vtu_grid)%cells%type_(i_node)
            if (vtu_grids(i_vtu_grid)%cells%type_(i_node)/=10) then
               write(uerr,*) 'Error in the program unit ',                     &
                  '"vtu_variable_reading". The ".vtu" cell type of the ',      &
                  'cell ',i_node,' is ',                                       &
                  vtu_grids(i_vtu_grid)%cells%type_(i_node),' instead of 10.'
            endif
         case ("connectivity")
            read(token,*)                                                      &
               vtu_grids(i_vtu_grid)%cells%points_IDs(i_node,i_component)
         case default
! Already treated in the calling program unit
      endselect
   enddo
enddo
!------------------------
! Deallocations
!------------------------
end subroutine vtu_variable_reading
#endif
