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
! Program unit: vtu_file_reading
! Description: Decoding of a ".vtu" file. Only tetrahedral cells are admitted 
!              (".vtu" cell type: 10).
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Further Copyright acknowledgments
! This first version of this program unit represented a modification of an 
! extract of the software tool vtu2FEHM (RSE SpA). Please refer to the git log 
! file for further details.
!-------------------------------------------------------------------------------
#if (defined SPACE_3D) && (defined SOLID_BODIES)
subroutine vtu_file_reading(i_vtu_grid,file_name,n_vtu_cells)
!------------------------
! Modules
!------------------------
use Hybrid_allocation_module
use Memory_I_O_interface_module
use Static_allocation_module
use I_O_file_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: i_vtu_grid
character(100),intent(in) :: file_name
integer(4),intent(inout) :: n_vtu_cells
integer(4) :: n_vtu_points,n_nodes,n_components,n_var_read,i_cell,i_tok,io_err
integer(4) :: i_sym,n_values,n_lines_var,n_tokens_line,ii
character(100) :: array_name,token,var_name,aux_char
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
   subroutine vtu_variable_reading(i_vtu_grid,n_nodes,n_components,n_lines_var,&
      n_tokens_line,var_name)
      implicit none
      integer(4),intent(in) :: i_vtu_grid,n_nodes,n_components,n_lines_var
      integer(4),intent(in) :: n_tokens_line
      character(100),intent(in) :: var_name
   end subroutine vtu_variable_reading
end interface
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
n_nodes = 0
n_var_read = 0
!------------------------
! Statements
!------------------------
! Open the ".vtu" file
call open_close_file(.true.,uvtu,file_name)
! Read the numbers of points and cells
do_points_cells: do
! Search the line containing the numbers of points and cells
   call ReadRiga(uvtu,line,io_err)
   if (io_err>0) then
      write(uerr,*) 'Error in reading the file ',file_name,' in the program ', &
         'unit "vtu_file_reading". The execution stops here.'
      stop
      elseif (io_err<0) then
! End of file
         exit
   endif
   i_tok = 2
   token = GetToken(i_tok,line,io_err)
   if (io_err/=0) cycle
   if (token(1:16)=='NumberOfPoints="') then
! Read the number of points 
      aux_char = trim(adjustl(token(17:100)))
      do i_sym=1,13
         if (aux_char(i_sym:i_sym)=='"') aux_char(i_sym:i_sym) = " "
      enddo
      read(aux_char,*) n_vtu_points
! Search the number of cells (token per token) in the remaining part of the line
      do
         i_tok = i_tok + 1 
         token = GetToken(i_tok,line,io_err)
         if (io_err==0) then
            if (token(1:15)=='NumberOfCells="') then
! Read the number of cells
               aux_char = trim(adjustl(token(16:100)))
               do i_sym=1,13
                  if (aux_char(i_sym:i_sym)=='"') aux_char(i_sym:i_sym) = " "
               enddo
               read (aux_char,*) n_vtu_cells
               exit do_points_cells
            endif
            else
! Error in reading the number of cells
               exit do_points_cells
         endif
      enddo
   endif
enddo do_points_cells
if ((n_vtu_points==0).or.(n_vtu_cells==0)) then
   write(uerr,*) 'The ".vtu" file',file_name,'read by the program unit ',      &
      '"vtu_file_reading" has no cell or no point. The program stops here.'
   stop
endif
write(ulog,*) 'The ".vtu" mesh n. ',i_vtu_grid,' is composed by ', n_vtu_points,  &
   ' ".vtu" points and ',n_vtu_cells,' ".vtu" cells.'
! Allocation and initialization of the ".vtu" global arrays
array_name = "vtu_cells_type"
call allocate_de_int4_r1(.true.,vtu_grids(i_vtu_grid)%cells%type_,n_vtu_cells, &
   array_name,ulog_flag=.true.)
vtu_grids(i_vtu_grid)%cells%type_(1:n_vtu_cells) = 0
array_name = "vtu_cells_points_IDs"
call allocate_de_int4_r2(.true.,vtu_grids(i_vtu_grid)%cells%points_IDs,        &
   n_vtu_cells,4,array_name,ulog_flag=.true.)
vtu_grids(i_vtu_grid)%cells%points_IDs(1:n_vtu_cells,1:4) = 0
! Allocation and initialization of the cell variable "surface"
! Cell faces are initialized as surface faces: it is sufficient to detect an 
! inner point to make them inner faces.
array_name = "vtu_cells_surface"
call allocate_de_log_r2(.true.,vtu_grids(i_vtu_grid)%cells%surface,            &
   n_vtu_cells,4,array_name,ulog_flag=.true.)
vtu_grids(i_vtu_grid)%cells%surface(1:n_vtu_cells,1:4) = .true.
array_name = "vtu_points_vertices"
call allocate_de_vertex_r1(.true.,vtu_grids(i_vtu_grid)%points%vertex,         &
   n_vtu_points,array_name,ulog_flag=.true.)
do ii=1,3
   vtu_grids(i_vtu_grid)%points%vertex(1:n_vtu_points)%pos(ii) = 0.d0
enddo
array_name = "vtu_points_surface"
call allocate_de_log_r1(.true.,vtu_grids(i_vtu_grid)%points%surface,           &
   n_vtu_points,array_name,ulog_flag=.true.)
vtu_grids(i_vtu_grid)%points%surface(1:n_vtu_points) = .false.
! Read the remaining lines of the file until the end of the file or the end of 
! reading of the useful variables. Lines containing "/" are not relevant: no 
! issue due to the the fact that anything after "/" on the same line is 
! automatically skipped.
do
! Check on the end of reading of the useful variables
   if (n_var_read==4) exit 
   call ReadRiga(uvtu,line,io_err)
   if (io_err>0) then
      write(uerr,*) 'Error in reading the variables in the file ',file_name,   &
         ' by the program unit "vtu_file_reading". The execution stops here.'
      stop
      elseif (io_err<0) then
! End of file
         exit
   endif
! Reading the first token of each label of the ".vtu" file sections
   token = GetToken(1,line,io_err)
! Recall the number of nodes of the current variable, or read its name and then 
! read/skip its values
   select case(trim(adjustl(token)))
      case("<PointData>","<Points>")
         n_nodes = n_vtu_points
      case("<CellData>","<Cells>")
         n_nodes = n_vtu_cells
      case("</PointData>","</Points>","</CellData>","</Cells>")
         n_nodes = 0
      case("<DataArray")
         i_tok = 2
         var_name = ""
         n_components = 1
! Loop over the following tokens of the same line
         do
            token = GetToken(i_tok,line,io_err)
            if (io_err/=0) exit
            if (token(1:6)=='Name="') then
! Read the variable name
               aux_char = trim(adjustl(token(7:100)))
               do i_sym=1,100
                  if (aux_char(i_sym:i_sym)=='"') aux_char(i_sym:i_sym) = " "
               enddo
               read(aux_char,'(a)') var_name
            endif
            if (token(1:20)=='NumberOfComponents="') then
! Read the variable components, if declared
               aux_char = trim(adjustl(token(21:100)))
               aux_char(2:2) = " "
               read(aux_char,*) n_components
            endif
            i_tok = i_tok + 1
         enddo
! Reading error
         if (trim(adjustl(var_name))=="") then
            write(uerr,*) 'The ".vtu" file n. ',i_vtu_grid,' read by the ',    &
               'program unit "vtu_file_reading" has found no name for the ',   &
               'current variable. The program stops here.'
            stop
         endif
! Read the first line of the values of the current variable
         call ReadRiga(uvtu,line,io_err)
! Assessment of the number of tokens in the current line
         i_tok = 0
         do while (io_err==0)
            i_tok = i_tok + 1
            token = GetToken(i_tok,line,io_err)
         enddo
         n_tokens_line = i_tok - 1
         if (n_tokens_line==0) then
            write(uerr,*) 'Program unit "vtu_file_reading". There is no value',&
               'for the variable ',trim(adjustl(var_name)),' .'
         endif
! Assessment of the number of lines of data of the current variable 
! depending on "n_nodes", "n_components" and "n_tokens_line"
         if (trim(adjustl(var_name))=="connectivity") n_components = 4
         n_values = n_nodes * n_components
         n_lines_var = int(n_values / n_tokens_line) + 1
         if (mod(n_values,n_tokens_line)==0) n_lines_var = n_lines_var - 1
! Step back one line
         backspace(uvtu)
! Read or skip the variable values depending on the variable name
         select case(trim(adjustl(var_name)))
            case("Points","nSurfaceLayers")
! The variable "nSurfaceLayers" can also be a cell-centred attribute
               if (n_nodes==n_vtu_points) then
                  call vtu_variable_reading(i_vtu_grid,n_nodes,n_components,   &
                     n_lines_var,n_tokens_line,var_name)
                  n_var_read = n_var_read + 1
                  else
                      write(aux_char,*) n_lines_var
                      aux_char = "(" // trim(adjustl(aux_char)) // "/)"
                      read(uvtu,aux_char)
               endif
            case("types")
               call vtu_variable_reading(i_vtu_grid,n_nodes,n_components,      &
                  n_lines_var,n_tokens_line,var_name)
               n_var_read = n_var_read + 1
            case("connectivity")
               call vtu_variable_reading(i_vtu_grid,n_nodes,n_components,      &
                  n_lines_var,n_tokens_line,var_name)
               n_var_read = n_var_read + 1
            case default
! Useless variable: the reading is skipped
               write(aux_char,*) n_lines_var
               aux_char = "(" // trim(adjustl(aux_char)) // "/)"
               read(uvtu,aux_char)
         endselect
      case default
! Keep reading in case other types of lines are encountered
   endselect
enddo
! The first ".vtu" vertex ID is nought: revision of the vertex IDs of each cell
!$omp parallel do default(none)                                                &
!$omp shared(n_vtu_cells,vtu_grids,i_vtu_grid)                                 &
!$omp private(i_cell)
do i_cell=1,n_vtu_cells
   vtu_grids(i_vtu_grid)%cells%points_IDs(i_cell,1:4) =                        &
      vtu_grids(i_vtu_grid)%cells%points_IDs(i_cell,1:4) + 1
enddo
!$omp end parallel do
! Close the ".vtu" file
call open_close_file(.false.,uvtu,file_name)
!------------------------
! Deallocations
!------------------------
! Deallocation of the array of the ".vtu" cell types
array_name = "vtu_cells_type"
call allocate_de_int4_r1(.false.,vtu_grids(i_vtu_grid)%cells%type_,            &
   array_name=array_name,ulog_flag=.true.)
end subroutine vtu_file_reading
#endif
