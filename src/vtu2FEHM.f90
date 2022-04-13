!-------------------------------------------------------------------------------
! "vtu2FEHM v.1.0" 
! Copyright 2016,2018 (RSE SpA)
! "vtu2FEHM v.1.0" authors and email contact are provided on 
! the documentation file.
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Description. “vtu2FEHM v.1.0” (RSE SpA) reads an unstructured grid (.vtu) of 
!              a finite element mesh, converts it into a grid input file for 
!              FEHM (.grid), produces the associated geological zone list 
!              (portion of the file .dat) and writes the node-elements list. 
!              Note: the mesh elements should be all of the same type 
!              (element%kind). Algorithm:
!                 1) Decoding the unstructured grid (.vtu)
!                 2) Writing the grid input file for FEHM (.grid)
!                 3) Writing the geological zone list (portion of the file .dat)
!                 4) Writing the unstructured .vtu file as a check
!                 5) Writing the node elements (starting from the element nodes)
!-------------------------------------------------------------------------------
PROGRAM vtu2FEHM
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
type node_type
   real :: coor(3)
   integer(4) :: material
   integer(4) :: element_ID(8)
end type node_type
type element_type
   integer(4) :: offset,kind,id,material
   integer(4) :: node_id(8)   
end type element_type
integer(4) :: n_materials,i,j,k,n_nodes,n_elements,element_nodes,ia,ia_pre
integer(4) :: alloc_stat,dealloc_stat,open_stat,close_stat
character(len=80) :: line,char1,char2,vtu_file_name,node_elements_file_name
integer(4),dimension(:),allocatable :: material_elem,material_nodes
character(len=4) :: fmt_char1,fmt_char2
type(node_type),dimension(:),allocatable :: node
type(element_type),dimension(:),allocatable :: element
character(len=80) :: vtu_file_name_test
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
write(*,*) " Start program Mvtu2FEHM"  
call getarg(1,vtu_file_name)
ia_pre=0
n_materials=0
char1(1:80)=" "
char2(1:80)=" "
vtu_file_name_test = "test.vtu"
node_elements_file_name = "node_elements.txt"
!------------------------
! Statements
!------------------------
! 1) Start decoding the unstructured grid (.vtu)
write(*,*) " 1) Decoding the unstructured grid (.vtu): start"  
open(11,file=vtu_file_name,IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) 'Opening of ',vtu_file_name,' failed; the program stops here. '
   stop
   else
      write(*,*) 'Opening of ',vtu_file_name,' successfully completed.'
endif
! General parameters
read(11,'(2/)')
read(11,'(a)') line
j=1
k=0
do i=1,len(line)
   ia=iachar(line(i:i))
   if ((ia>=48).and.(ia<=57)) then
      if ((ia_pre>=48).and.(ia_pre<=57)) then
         j=j+1
         else
            k=k+1
            j=1
      endif
      if (k==1) then
         char1(j:j)=achar(ia)
         else   
            char2(j:j)=achar(ia)
      endif
   endif
   ia_pre=ia
enddo
fmt_char1='(i )'
fmt_char2='(i )'
write(fmt_char1(3:3),'(i1)') len_trim(char1)
write(fmt_char2(3:3),'(i1)') len_trim(char2)
read(char1,fmt_char1) n_nodes 
read(char2,fmt_char2) n_elements 
if(.not.allocated(node)) then
   allocate(node(n_nodes),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) 'Allocation of "node" failed; the program stops here. '
      stop
      else
         write(*,*) 'Allocation of "node" successfully completed.'
   endif
endif
if(.not.allocated(element)) then
   allocate(element(n_elements),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) 'Allocation of "element" failed; the program stops here. '
      stop
      else
         write(*,*) 'Allocation of "element" successfully completed.'
   endif
endif
write(*,*) "        The mesh is composed by ", n_nodes, "nodes and ",          &
   n_elements," elements"
! Reading the nodes 
read(11,'(/)')
do i=1,n_nodes,12
   read(11,'(1x,12(3(f12.4),2x))') (node(j)%coor(:),j=i,i+11)
enddo
read(11,'(2/)')
! Reading the connectivity structure        
read(11,'()')
do i=1,(n_elements-7),7
   read(11,'(56(1x,i7))') (element(j)%node_id(:),j=i,i+6)
enddo
do k=i,n_elements
   read(11,'(8(1x,i7))',advance='no') element(k)%node_id(:)
enddo
read(11,'(/)') 
! Reading the offsets (position of the last element node in the .vtu list)
read(11,'(a)')   
do i=1,(n_elements-56),56
   read(11,'(56(1x,i7))') (element(j)%offset,j=i,i+55)
enddo
do k=i,n_elements
   read(11,'(1x,i7)',advance='no') element(k)%offset
enddo
read(11,'(/)')
! Reading the element type
read(11,'()')    
do i=1,(n_elements-96),96
   read(11,'(96(1x,i2))') (element(j)%kind,j=i,i+95)
enddo
do k=i,n_elements
   read(11,'(1x,i2)',advance='no') element(k)%kind
enddo
read(11,'(3/)')
! Reading the mesh element indeces
read(11,'()')
do i=1,(n_elements-56),56
   read(11,'(56(1x,i7))') (element(j)%id,j=i,i+55)
enddo
do k=i,n_elements
   read(11,'(1x,i7)',advance='no') element(k)%id
enddo
read(11,'(/)')
! Reading the materials
read(11,'()')
do i=1,(n_elements-96),96
   read(11,'(96(1x,i2))') (element(j)%material,j=i,i+95)
enddo
do k=i,n_elements
   read(11,'(1x,i2)',advance='no') element(k)%material
enddo
do i=1,n_elements
   n_materials=max(n_materials,element(i)%material)
enddo
if(.not.allocated(material_elem)) then
   allocate(material_elem(n_materials),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) 'Allocation of "material_elem" failed; the program stops here.'
      stop
      else
         write(*,*) 'Allocation of "material_elem" successfully completed.'
   endif
endif
do i=1,n_materials
   material_elem(i)=0
enddo
! Nota: +1 perhè il primo ID dei vtu è zero (convenzione Phyton)
do i=1,n_elements
   do j=1,8
      element(i)%node_id(j)=element(i)%node_id(j)+1
   enddo
enddo
flush(11)
close(11,IOSTAT=close_stat)
if (close_stat/=0) then
   write(*,*) 'Closing of ',vtu_file_name,' failed; the program stops here. '
   stop
   else
      write(*,*) 'Closing of ',vtu_file_name,' successfully completed.'
endif
write(*,*) "       Decoding the unstructured grid (.vtu): end" 
! End decoding the unstructured grid (.vtu)
! 2) Writing the grid input file for FEHM (.grid)
write(*,*) " 2) Writing the grid input file for FEHM (.grid): start"  
open(12,file="FEHM_mesh.grid",IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) 'Opening of FEHM_mesh.grid failed; the program stops here. '
   stop
   else
      write(*,*) 'Opening of FEHM_mesh.grid successfully completed.'
endif
write(12,'(a)') "coor   n/a     "
write(12,'(i12)') n_nodes
j=1
do i=1,n_nodes
write(12,'(i10,3(f15.5))') i,node(i)%coor(:)
enddo
write(12,'(a)') "         0        0.00000        0.00000        0.00000" 
select case (element(1)%kind)
   case (12)
   element_nodes=8
end select
write(12,'(a)') "elem trad"
write(12,'(2(i8))',advance='NO') element_nodes,n_elements
write(12,'(a)') "   0"
do i=1,n_elements
   write(12,'(9(i8))') i,element(i)%node_id(5),element(i)%node_id(6),          &
                       element(i)%node_id(7),element(i)%node_id(8),            &
                       element(i)%node_id(1),element(i)%node_id(2),            &
                       element(i)%node_id(3),element(i)%node_id(4)
   material_elem(element(i)%material)=material_elem(element(i)%material)+1
enddo
write(12,'(a)') "       0       0       0       0       0       0       0       0       0"
write(12,'(a)') "stop"
flush(12)
close(12,IOSTAT=close_stat)
if (close_stat/=0) then
   write(*,*) 'Closing of FEHM_mesh.grid failed; the program stops here. '
   stop
   else
      write(*,*) 'Closing of FEHM_mesh.grid successfully completed.'
endif
write(*,*) "       Writing the grid input file for FEHM (.grid): end" 
! Writing the grid input file for FEHM (.grid): end
! 3) Writing the geological zone list (portion of the file .dat)
write(*,*) " 3) Writing the geological zone list (portion of the file .dat",   &
   "): start" 
! Each node needs to be given a material
do i=1,n_nodes
   node(i)%material = 0
   do j=1,n_elements
      do k=1,8 
         if (element(j)%node_id(k)==i) then
            node(i)%material = max(node(i)%material,element(j)%material)
         endif  
      enddo
   enddo
enddo
if(.not.allocated(material_nodes)) then
   allocate(material_nodes(n_materials),STAT=alloc_stat)
   if (alloc_stat/=0) then
      write(*,*) 'Allocation of "material_nodes" failed; the program stops ',  &
         'here.'
      stop
      else
         write(*,*) 'Allocation of "material_nodes" successfully completed.'
   endif
endif
do i=1,n_materials
   material_nodes(i)=0
enddo
do i=1,n_nodes
   material_nodes(node(i)%material)=material_nodes(node(i)%material)+1
enddo
open(13,file="FEHM_zones.dat",IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) 'Opening of FEHM_zones.dat failed; the program stops here. '
   stop
   else
      write(*,*) 'Opening of FEHM_zones.dat successfully completed.'
endif
do i=1,n_materials
   write(13,'(i8)') i
   write(13,'(a)') "nnum"
   write(13,'(i8)') material_nodes(i)
   k=1
   do j=1,n_nodes
      if (node(j)%material==i) then
         k=k+1
         if (mod(k,56)/=0) then
            write(13,'(i8)',ADVANCE='NO') j
            else
               write(13,'(i8)') j
         endif
      endif
   enddo
   write(13,'()')
enddo
flush(13)
close(13,IOSTAT=close_stat)
if (close_stat/=0) then
   write(*,*) 'Closing of FEHM_zones.dat failed; the program stops here. '
   stop
   else
      write(*,*) 'Closing of FEHM_zones.dat successfully completed.'
endif
write(*,*) "       Writing the geological zone list (portion of the file .dat",&
   "): end"
! Writing the geological zone list (portion of the file .dat): end
! 4) Writing the unstructured .vtu file as a check 
write(*,*) " 4) Writing the unstructured .vtu file as a check: start  "  
open(14,file=vtu_file_name_test,IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) 'Opening of ',vtu_file_name_test,' failed; the program stops',   &
      'here. '
   stop
   else
      write(*,*) 'Opening of ',vtu_file_name_test,' successfully completed.'
endif
! General parameters
write(14,'(a)') '<?xml version="1.0"?>'
write(14,'(a)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
write(14,'(a)') '<UnstructuredGrid>'
write(14,'(a)') line
write(14,'(a)') '    <Points>'
write(14,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
! writing the nodes 
do i=1,n_nodes,12
   write(14,'(1x,12(3(f12.4),2x))') (node(j)%coor(:),j=i,i+11)
enddo
write(14,'(a)') '      </DataArray>'
! Writing the connectivity structure        
write(14,'(a)') '    </Points>'
write(14,'(a)') '    <Cells>'
write(14,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
do i=1,(n_elements-7),7
   write(14,'(56(1x,i7))') (element(j)%node_id(:)-1,j=i,i+6)
enddo
do k=i,n_elements
   write(14,'(8(1x,i7))',advance='no') element(k)%node_id(:)-1
enddo
write(14,'()')
write(14,'(a)') '      </DataArray>'  
! Writing the offsets (position of the last element node in the .vtu list)
write(14,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
do i=1,(n_elements-56),56
   write(14,'(56(1x,i7))') (element(j)%offset,j=i,i+55)
enddo
do k=i,n_elements
   write(14,'(1x,i7)',advance='no') element(k)%offset
enddo
write(14,'()')
write(14,'(a)') '      </DataArray>'
! Writing the element type
write(14,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
do i=1,(n_elements-96),96
   write(14,'(96(1x,i2))') (element(j)%kind,j=i,i+95)
enddo
do k=i,n_elements
   write(14,'(1x,i2)',advance='no') element(k)%kind
enddo
write(14,'()')
write(14,'(a)') '      </DataArray>'
! Writing the mesh element indeces
write(14,'(a)') '    </Cells>'
write(14,'(a)') '    <CellData>'
write(14,'(a)') '      <DataArray type="Int32" Name="MeshIndex" format="ascii" >'
do i=1,(n_elements-56),56
   write(14,'(56(1x,i7))') (element(j)%id,j=i,i+55)
enddo
do k=i,n_elements
   write(14,'(1x,i7)',advance='no') element(k)%id
enddo
write(14,'()')
write(14,'(a)') '      </DataArray>'
! Writing the materials
write(14,'(a)') '      <DataArray type="Int32" Name="Materials" format="ascii" >'
do i=1,(n_elements-96),96
   write(14,'(96(1x,i2))') (element(j)%material,j=i,i+95)
enddo
do k=i,n_elements
   write(14,'(1x,i2)',advance='no') element(k)%material
enddo
write(14,'()')
write(14,'(a)') '      </DataArray>'
write(14,'(a)') '    </CellData>'
write(14,'(a)') '    <PointData>'
write(14,'(a)') '    </PointData>'
write(14,'(a)') '  </Piece>'
write(14,'(a)') '</UnstructuredGrid>'
write(14,'(a)') '</VTKFile>'
flush(14)
close(14,IOSTAT=close_stat)
if (close_stat/=0) then
   write(*,*) 'Closing of ',vtu_file_name_test,' failed; the program stops ',  &
      'here. '
   stop
   else
      write(*,*) 'Closing of ',vtu_file_name_test,' successfully completed.'
endif
write(*,*) "       Writing the unstructured .vtu file as a check : end"
! Writing the unstructured .vtu file as a check: end 
! 5) Writing the node elements (starting from the element nodes): start
write(*,*) " 5) Writing the node elements (starting from the element nodes): ",&
   "start  "
open(15,file=node_elements_file_name,IOSTAT=open_stat)
if (open_stat/=0) then
   write(*,*) 'Opening of ',node_elements_file_name,' failed; the program ',   &
      'stops here. '
   stop
   else
      write(*,*) 'Opening of ',node_elements_file_name,' successfully ',       &
         'completed.'
endif
do i=1,n_nodes
   node(i)%element_ID(:) = -999
enddo
! Loop over the elements
do i=1,n_elements
! Loop over the element nodes
   do j=1,8
! Loop over the node elements
      do k=1,8
         if (node(element(i)%node_id(j))%element_ID(k)>=0) then
            cycle
            else
               node(element(i)%node_id(j))%element_ID(k) = i - 1
               exit
         endif
      enddo
   enddo
enddo
write(15,'(a)') "Node ID Element IDs"
do i=1,n_nodes
   write(15,'(9(1x,i7))') (i-1),node(i)%element_ID(:)
enddo
flush(15)
close(15,IOSTAT=close_stat)
if (close_stat/=0) then
   write(*,*) 'Closing of ',node_elements_file_name,' failed; the program ',   &
      'stops here. '
   stop
   else
      write(*,*) 'Closing of ',node_elements_file_name,' successfully ',       &
         'completed.'
endif
write(*,*) "       Writing the node elements (starting from the element nodes",&
   ") : end"
! Writing the node elements (starting from the element nodes): end
!------------------------
! Deallocations
!------------------------
write(*,*) " End program vtu2FEHM"
end program vtu2FEHM

