!AA601 the whole subroutine
!cfile BC_wall_elements.f90
!************************************************************************************
!                          S P H E R A 5.0.3 
!
!                Smoothed Particle Hydrodinamics Code
!
!************************************************************************************
!
! Subroutine : DBSPH_IC_surface_elements
!
! Creation   : Amicarelli A., 16Jan15; Main features (see Purpose)
!
!************************************************************************************
! Purpose            : Initialization of wall surface elements for DBSPH
! Calling procedures : gest_input
! Called procedures  : /
!************************************************************************************

subroutine DBSPH_IC_surface_elements

! Modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
use DIAGNOSTIC_MODULE

! Declarations
implicit none
integer(4) :: alloc_stat,i,j
integer(4),external :: ParticleCellNumber

! Interface blocks

! Allocations

! Initializations

! Statements
if (.not.allocated(pg_w)) then
   allocate(pg_w(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet),STAT=alloc_stat) 
   if (alloc_stat/=0) then
      write(nout,*) 'Allocation of pg_w in DBSPH_IC_surface_elements failed; the program terminates here.'
      call diagnostic (arg1=5,arg2=340)
      stop ! Stop the main program
      else
         pg_w(:)%cella = 0
         pg_w(:)%adjacent_faces(1) = 0
         pg_w(:)%adjacent_faces(2) = 0
         pg_w(:)%adjacent_faces(3) = 0         
         pg_w(:)%coord(1) = 0.d0
         pg_w(:)%coord(2) = 0.d0
         pg_w(:)%coord(3) = 0.d0
         pg_w(:)%vel(1) = 0.d0
         pg_w(:)%vel(2) = 0.d0         
         pg_w(:)%vel(3) = 0.d0
         pg_w(:)%dens = 0.d0
         pg_w(:)%pres = 0.d0
         pg_w(:)%normal(1) = 0.d0
         pg_w(:)%normal(2) = 0.d0
         pg_w(:)%normal(3) = 0.d0
         pg_w(:)%weight = 0.d0 
         pg_w(:)%wet = 0 
         pg_w(:)%mass = 0.d0 
         pg_w(:)%k_d = 0.d0
         pg_w(:)%volume = 0.d0
         write (nout,*) "Allocation of pg_w in DBSPH_IC_surface_elements successfully completed."
   endif   
endif 
!$omp parallel do default(none) shared(DBSPH,pg_w,Med,pg,ncord) private(i)
do i=1,DBSPH%n_w 
!AA601!!! sub start
   if (ncord==3) then
      pg_w(i)%coord(:) = ( DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(1))%pos(:) + &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(2))%pos(:) + &
                           DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(3))%pos(:) ) / 3.d0
      else
         pg_w(i)%coord(:) = ( DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(1))%pos(:) + &
                              DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(2))%pos(:) + &
                              DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(3))%pos(:) + &
                              DBSPH%surf_mesh%vertices(DBSPH%surf_mesh%faces(i)%vert_list(4))%pos(:) ) / 4.d0    
   endif
   pg_w(i)%cella = ParticleCellNumber(pg_w(i)%coord)
   pg_w(i)%vel(:) = 0.d0
   pg_w(i)%dens = Med(1)%den0                         
   pg_w(i)%pres = 0.d0                          
   pg_w(i)%normal(:) = DBSPH%surf_mesh%faces(i)%normal(:)                      
   pg_w(i)%weight = DBSPH%surf_mesh%faces(i)%area
   pg_w(i)%wet = 0       
   pg_w(i)%adjacent_faces(:) = 0
end do
!$omp end parallel do
!Initializing fictitious surface elements representing DBSPH inlet sections
!$omp parallel do default(none) shared(DBSPH,pg_w,Med,pg) private(i,j)
do i=(DBSPH%n_w+1),(DBSPH%n_w+DBSPH%n_inlet)
   j= i-DBSPH%n_w  
   pg_w(i)%coord(:) = DBSPH%inlet_sections(j,1:3)
   pg_w(i)%cella = ParticleCellNumber(pg_w(i)%coord)
   pg_w(i)%vel(:) = DBSPH%inlet_sections(j,7:9)
   pg_w(i)%dens = Med(1)%den0                         
   pg_w(i)%pres = 0.d0                          
   pg_w(i)%normal(:) = DBSPH%inlet_sections(j,4:6)                      
   pg_w(i)%weight = 0.d0
   pg_w(i)%wet = 0       
   pg_w(i)%adjacent_faces(:) = 0
end do
!$omp end parallel do
!Initializing fictitious surface elements representing DBSPH outlet sections
!$omp parallel do default(none) shared(DBSPH,pg_w,Med,pg) private(i,j)
do i=(DBSPH%n_w+DBSPH%n_inlet+1),(DBSPH%n_w+DBSPH%n_inlet+DBSPH%n_outlet)
   j= i-(DBSPH%n_w+DBSPH%n_inlet)  
   pg_w(i)%coord(:) = DBSPH%outlet_sections(j,1:3)
   pg_w(i)%cella = ParticleCellNumber(pg_w(i)%coord)
   pg_w(i)%vel(:) = 0.d0
   pg_w(i)%dens = Med(1)%den0                         
   pg_w(i)%pres = 0.d0                          
   pg_w(i)%normal(:) = DBSPH%outlet_sections(j,4:6)                      
   pg_w(i)%weight = 0.d0
   pg_w(i)%wet = 0       
   pg_w(i)%adjacent_faces(:) = 0
end do
!$omp end parallel do
! Deallocations

return
end subroutine DBSPH_IC_surface_elements
!end of the subroutine

