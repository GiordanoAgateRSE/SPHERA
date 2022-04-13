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
! Program unit: vtu_reading
! Description: Decoding of a ".vtu" file
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Further Copyright acknowledgments
! This program unit represents a modification of an extract of the software 
! tool vtu2FEHM (RSE SpA). Please refer to the git log file for further 
! details.
!-------------------------------------------------------------------------------
subroutine vtu_reading
!------------------------
! Modules
!------------------------
!------------------------
! Declarations
!------------------------
implicit none
type node_type
   real :: coor(3)
end type node_type
type element_type
   integer(4) :: offset,ekind,id,material
   integer(4) :: node_id(8)   
end type element_type
integer(4) :: i,j,k,n_nodes,n_elements,ia,ia_pre
integer(4) :: alloc_stat,dealloc_stat,open_stat,close_stat
character(len=80) :: line,char1,char2,vtu_file_name
character(len=4) :: fmt_char1,fmt_char2
type(node_type),dimension(:),allocatable :: node
type(element_type),dimension(:),allocatable :: element
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
ia_pre=0
char1(1:80)=" "
char2(1:80)=" "
!------------------------
! Statements
!------------------------
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
   read(11,'(96(1x,i2))') (element(j)%ekind,j=i,i+95)
enddo
do k=i,n_elements
   read(11,'(1x,i2)',advance='no') element(k)%ekind
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
! The first ID is nought
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
!------------------------
! Deallocations
!------------------------
end subroutine vtu_reading
