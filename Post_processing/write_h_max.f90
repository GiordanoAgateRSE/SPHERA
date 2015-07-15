!cfile write_h_max.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : write_h_max
!
! Versions: 
! 01 Amicarelli 08Apr14 (creation)
!
!************************************************************************************
! Module purpose : (v5.04) Compute and write the 2D array of the maximum values of the water
!                  depth, at the nodes of a Cartesian topography (only in 3D),
!                  together with the 2D field of the maximum (over time) specific 
!                  flow rates.
!
! Calling routine: Gest_Trans
!
! Called routines: /
!
!************************************************************************************

subroutine write_h_max

!Using modules
 use GLOBAL_MODULE
 use AdM_USER_TYPE
 use ALLOC_MODULE
 use files_entities
 
!Declarations
 implicit none
 integer(4) :: i_zone,i_vertex,GridColumn
 double precision :: pos(3)
 double precision,dimension(:),allocatable :: h_max 
 character(255) :: nomefile_h_max
 
!External function
 integer(4),external :: ParticleCellNumber
 
! h_max .txt file: creation and heading 
 write(nomefile_h_max,"(a,a)") nomecaso(1:len_trim(nomecaso)),"_h_max_q_max.txt"
 open(ncpt,file=nomefile_h_max,status="unknown",form="formatted")
 write(ncpt,*) "Maximum water depth (m) and specific flow rate (m^2/s)"
 write(ncpt,'(6(a))') "           x(m)","           y(m)","       h_max(m)"," Z_fluid_max(m)","     z_topog(m)","   q_max(m^2/s)"
 flush(ncpt) 
 
!Statements
 do i_zone=1,NPartZone
    if (Partz(i_zone)%IC_source_type == 2) then
!Allocating h_max 
       allocate(h_max(Partz(i_zone)%npoints))
!Initializing h_max
       h_max = 0.
!$omp parallel do default(none) shared(Partz,Vertice,Grid,h_max,Z_fluid_max,ncpt,i_zone,q_max) private(i_vertex,GridColumn,pos)
       do i_vertex=Partz(i_zone)%ID_first_vertex,Partz(i_zone)%ID_last_vertex
          pos(1) = Vertice(1,i_vertex)
          pos(2) = Vertice(2,i_vertex)
          pos(3) = Grid%extr(3,1)+0.0000001
          GridColumn = ParticleCellNumber(pos)
          h_max(i_vertex-Partz(i_zone)%ID_first_vertex+1) = max ((Z_fluid_max(GridColumn)-Vertice(3,i_vertex)),0.d0)
!$omp critical (omp_write_h_max)
          write(ncpt,'(6(f14.4,1x))')Vertice(1,i_vertex),Vertice(2,i_vertex),h_max(i_vertex-Partz(i_zone)%ID_first_vertex+1),Z_fluid_max(GridColumn),Vertice(3,i_vertex) &
                                     ,q_max(i_vertex-Partz(i_zone)%ID_first_vertex+1)        
!$omp end critical (omp_write_h_max)
       end do
!$omp end parallel do    
       exit
    endif
 end do

! h_max .txt file: closing
 close (ncpt)

!Deallocations
 if (allocated(h_max)) deallocate(h_max)
 if (allocated(q_max)) deallocate(q_max)

 return
end subroutine write_h_max
!---split

