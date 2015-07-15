!cfile Update_Zmax_at_grid_vert_columns.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : Update_Zmax_at_grid_vert_columns
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Updating the 2D array of the maximum values of the fluid particle
!                  height, at the nodes of the grid vertical columns (only in 3D).
!                  Printing the 2D field of the water depth (current time step), 
!                  according to output frequency chosen in input file (only in 3D).
!                  Printing the 2D fields of specific flow rate components (current
!                  time step), at the same frequency of the water depth (only in 3D).
!
! Calling routine: Loop_Irre_3D
!
! Called routines: /
!
!************************************************************************************

subroutine Update_Zmax_at_grid_vert_columns(print_flag)

!Using modules
 use GLOBAL_MODULE
 use AdM_USER_TYPE
 use ALLOC_MODULE
 use files_entities

! Declarations
 implicit none
 integer(4) :: npi,GridColumn
 double precision :: pos(3)
 double precision, allocatable, dimension(:) :: Z_fluid_step,h_step,qx_step,qy_step,qx_step_grid,qy_step_grid,n_part_step 
 integer(4) :: i_zone,i_vertex,i_aux,i_grid,j_grid
 integer(4), intent(in) :: print_flag
 character(255) :: nomefile_h_step
 
!External function 
 integer(4),external :: ParticleCellNumber

!Allocations
 do i_zone=1,NPartZone
    if (Partz(i_zone)%IC_source_type == 2) then
       allocate(h_step(Partz(i_zone)%npoints))
       allocate(qx_step(Partz(i_zone)%npoints))
       allocate(qy_step(Partz(i_zone)%npoints))
       exit
    endif
 enddo
 allocate(Z_fluid_step(Grid%ncd(1)*Grid%ncd(2)))
 allocate(qx_step_grid(Grid%ncd(1)*Grid%ncd(2)))
 allocate(qy_step_grid(Grid%ncd(1)*Grid%ncd(2)))
 allocate(n_part_step(Grid%ncd(1)*Grid%ncd(2)))
 Z_fluid_step = -999.d0
 h_step = 0.d0
 qx_step = 0.d0
 qy_step = 0.d0
 qx_step_grid = 0.d0
 qy_step_grid = 0.d0 
 n_part_step = 0
 
! Statements
!$omp parallel do default(none) shared(nag,pg,Grid,Z_fluid_max,Domain,Z_fluid_step,n_part_step,qx_step_grid,qy_step_grid,tempo) private(npi,GridColumn,pos)
 do npi=1,nag
    pos(1) = pg(npi)%coord(1)
    pos(2) = pg(npi)%coord(2)
    pos(3) = Grid%extr(3,1) + 0.0000001
    GridColumn = ParticleCellNumber(pos)
    Z_fluid_step(GridColumn) = max(Z_fluid_step(GridColumn),(pg(npi)%coord(3)+0.5*Domain%dd)) 
    qx_step_grid(GridColumn) = qx_step_grid(GridColumn) + pg(npi)%vel(1)
    qy_step_grid(GridColumn) = qy_step_grid(GridColumn) + pg(npi)%vel(2)
    n_part_step(GridColumn) = n_part_step(GridColumn) + 1
    Z_fluid_max(GridColumn) = max (Z_fluid_max(GridColumn),(pg(npi)%coord(3)+0.5*Domain%dd))
 end do
!$omp end parallel do
!$omp parallel do default(none) shared(Grid,n_part_step,qx_step_grid,qy_step_grid,tempo) private(GridColumn,pos,i_grid,j_grid)
    do i_grid=1,Grid%ncd(1)
       do j_grid=1,Grid%ncd(2)
          pos(1) = (i_grid-0.5) * Grid%dcd(1) + Grid%extr(1,1)
          pos(2) = (j_grid-0.5) * Grid%dcd(2) + Grid%extr(2,1)
          pos(3) = Grid%extr(3,1) + 0.0000001
          GridColumn = ParticleCellNumber(pos)
          if (n_part_step(GridColumn)>0) then
             qx_step_grid(GridColumn) = qx_step_grid(GridColumn)/n_part_step(GridColumn)
             qy_step_grid(GridColumn) = qy_step_grid(GridColumn)/n_part_step(GridColumn)
             else
                qx_step_grid(GridColumn) = 0.d0
                qy_step_grid(GridColumn) = 0.d0
          endif
       enddo
    enddo
!$omp end parallel do
! .txt file creation and heading (only at the beginning of the simulation)
 if (it_corrente == 1) then
    write(nomefile_h_step,"(a,a)") nomecaso(1:len_trim(nomecaso)),"_h_qx_qy_step.txt"
    open(ncpt,file=nomefile_h_step,status="unknown",form="formatted")
    write(ncpt,*) "Water depth (h), Specific flow rate components (q_x,q_y)"
    write(ncpt,'(7(a))') "           x(m)","           y(m)","  h_step(m)"," Z_fl_max_stp(m)","     z_topog(m)","     q_x(m^2/s)","     q_y(m^2/s)"
    flush(ncpt)
    else 
! Writing the 2D free surface field at the current time step
       if (print_flag==1) then
          write(nomefile_h_step,"(a,a,i8.8,a)") nomecaso(1:len_trim(nomecaso)),'_h_qx_qy_step',it_corrente,".txt"
          open (ncpt,file=nomefile_h_step,status="unknown",form="formatted")
       endif   
       do i_zone=1,NPartZone
          if (Partz(i_zone)%IC_source_type == 2) then
!$omp parallel do default(none) shared(ncpt,Partz,Vertice,Grid,h_step,Z_fluid_step,i_zone,qx_step,qy_step,qx_step_grid,qy_step_grid,n_part_step,q_max,tempo,print_flag) &
!$omp private(i_vertex,GridColumn,pos,i_aux)
             do i_vertex=Partz(i_zone)%ID_first_vertex,Partz(i_zone)%ID_last_vertex
                i_aux = i_vertex-Partz(i_zone)%ID_first_vertex+1 
                pos(1) =  Vertice(1,i_vertex)
                pos(2) = Vertice(2,i_vertex)
                pos(3) = Grid%extr(3,1)+0.0000001
                GridColumn = ParticleCellNumber(pos)
                h_step(i_aux) = max ((Z_fluid_step(GridColumn)-Vertice(3,i_vertex)),0.d0)
                qx_step(i_aux) = qx_step_grid(GridColumn)
                qy_step(i_aux) = qy_step_grid(GridColumn)
                if (qx_step_grid(GridColumn)/=-999.d0) then
                    qx_step(i_aux) = qx_step_grid(GridColumn)*h_step(i_aux)
                    qy_step(i_aux) = qy_step_grid(GridColumn)*h_step(i_aux)
                    q_max(i_aux) = max(q_max(i_aux),dsqrt(qx_step(i_aux)**2 + qy_step(i_aux)**2))
                endif    
                if (print_flag==1) then
!$omp critical (omp_write_h_step)
                   write (ncpt,'(7(f14.4,1x))') Vertice(1,i_vertex),Vertice(2,i_vertex),h_step(i_aux),Z_fluid_step(GridColumn),Vertice(3,i_vertex), &
                                                qx_step(i_aux),qy_step(i_aux)
!$omp end critical (omp_write_h_step)
                endif
             end do
!$omp end parallel do    
             exit
          endif
       enddo
 endif

! Closing the step file
 if (print_flag==1) close (ncpt)
 
!Deallocations
 if (allocated(Z_fluid_step)) deallocate(Z_fluid_step)
 if (allocated(h_step)) deallocate(h_step)
 if (allocated(qx_step)) deallocate(qx_step)
 if (allocated(qy_step)) deallocate(qy_step)
 if (allocated(qx_step_grid)) deallocate(qx_step_grid)
 if (allocated(qy_step_grid)) deallocate(qy_step_grid)
 if (allocated(n_part_step)) deallocate(n_part_step)
 
 return
end subroutine Update_Zmax_at_grid_vert_columns
!---split

