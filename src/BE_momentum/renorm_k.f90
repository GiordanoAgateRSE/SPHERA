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
! Program unit: renorm_k
! Description:  Reormalization algorithm on the Gallati Anti Cluster Kernel 
! derivative.
!-------------------------------------------------------------------------------
subroutine renorm_k(npi,Binv)
!------------------------
! Modules
!------------------------
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(in) :: npi
integer(4) :: npj,contj,npartint,ii,index_rij_su_h
double precision, dimension (3,1) :: gradk
double precision, dimension (1,3) :: gradk_t
double precision, dimension(3,3) :: B,gradk_incr
double precision, dimension(2,2) :: B2
double precision, intent(out) ,dimension(3,3) :: Binv
double precision :: detB
double precision :: rijtemp,rij_su_h
double precision :: ragtemp(3)
#ifdef SPACE_3D
integer(4) :: test_m
#endif
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
gradk_t(:,:) = zero
gradk(:,:) = zero
gradk_incr(:,:) = zero
B(:,:) = zero
B2(:,:) = zero
Binv(:,:) = zero
detB = zero
!------------------------
! Statements
!------------------------
!loop over the neighbouring particles of i
do contj=1,nPartIntorno(npi)
   npartint = (npi - 1) * NMAXPARTJ + contj
   npj = PartIntorno(npartint)
! To check if this term is needed in the previous loop 
   if (npi==npj) cycle
   ! Relative distance
   ragtemp(1:3) = rag(:,npartint)
   ! Square of the relative distance (to improve accuracy later on)
                  rijtemp = ragtemp(1)*ragtemp(1) + ragtemp(2)*ragtemp(2) +       &
                            ragtemp(3)*ragtemp(3)
   rijtemp = dsqrt(rijtemp)
   rij_su_h = rijtemp / Domain%h  
   index_rij_su_h = int(rij_su_h)
   if (index_rij_su_h>=2) cycle
!trasposing vector column to row
    gradk(:,1) = -PartKernel(3,npartint)*rag(:,npartint)
    call MatrixTransposition(gradk(:,1),gradk_t(1,:),1,3)
!sum of the contributions over the neighbouring particles
    call MatrixProduct(rag(:,npartint),gradk_t(1,:),gradk_incr(:,:),3,1,3)
    B (:,:)= B(:,:)+ gradk_incr * (pg(npj)%mass/pg(npj)%dens)
enddo

#ifdef SPACE_2D
  B2(1,1) = B(1,1)
  B2(1,2) = B(1,3)
  B2(2,1) = B(3,1)
  B2(2,2) = B(3,3)
  detB = B2(1,1)*B2(2,2)-B2(1,2)*B2(2,1)
  if(abs(detB) >= input_any_t%Btol) then
! Calculate inverse of B2 (2x2 matrix)
    Binv(1,1)=B2(2,2)/detB
    Binv(3,1)=-B2(1,2)/detB
    Binv(1,3)=-B2(2,1)/detB
    Binv(3,3)=B2(1,1)/detB
  else
    do ii=1,3
      Binv(ii,ii)=1
    enddo
endif
#endif

#ifdef SPACE_3D
    test_m = 1
    detB = B(1,1) * B(2,2) * B(3,3) + B(2,1) * B(3,2) * B(1,3) +                  &
           B(3,1) * B(1,2) * B(2,3) - B(1,1) * B(3,2) * B(2,3) -                  &
           B(3,1) * B(2,2) * B(1,3) - B(2,1) * B(1,2) * B(3,3)    
if (abs(detB) >= input_any_t%Btol) then
    call Matrix_Inversion_3x3(B(:,:),Binv(:,:),test_m)
else
    do ii=1,3
        Binv(ii,ii)=1
    enddo
endif

#endif
!------------------------
! Deallocations
!------------------------
return
end subroutine renorm_k
