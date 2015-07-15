!AA501b the whole subroutine
!AA504 sub
!cfile reference_system_change.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! Module name     : reference_system_change
!
! Creation : Amicarelli-Agate, 24Oct12
!
!************************************************************************************
! Module purpose : Transformation of the coordinates, expressed in a new reference system
!
!AA601
! Calling routine: RHS_body_dynamics,DBSPH_inlet_outlet
!
! Called routines: MatrixTransposition 
!
!************************************************************************************

subroutine reference_system_change(pos_old_ref,new_origin_old_ref,new_cos_dir_old_ref,pos_new_ref)

! Declarations
 implicit none
 integer(4) :: i,aux_int
 double precision,intent(IN) :: pos_old_ref(3),new_origin_old_ref(3) 
 double precision,intent(IN) :: new_cos_dir_old_ref(3,3)
 double precision,intent(INOUT) :: pos_new_ref(3)
 double precision :: aux_vec(3)
 double precision :: inv_new_cos_dir_old_ref(3,3)
 
! Statements
 aux_vec(:) = pos_old_ref(:) - new_origin_old_ref(:)
 call Matrix_Inversion_3x3(new_cos_dir_old_ref,inv_new_cos_dir_old_ref,aux_int)
 do i=1,3
    pos_new_ref(i) = dot_product(inv_new_cos_dir_old_ref(i,:),aux_vec)
 enddo

return
end subroutine reference_system_change
!---split

