!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-Ricerca di Sistema)



! SPHERA authors and email contact are provided on SPHERA documentation.

! This file is part of SPHERA v.8.0.
! SPHERA v.8.0 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
! You should have received a copy of the GNU General Public License
! along with SPHERA. If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: ComputeBoundaryIntegralTab                                   
! Description: To compute local coordinates (x,y,z) of a grid of points, regularly distributed on the semisphere z<0 
!              (radius = 2h), whose centre is the origin O of local axis. The semisphere will be superposed to the influence 
!              sphere (kernel support) of a generic particle near a plane boundary face, and oriented in such a way that
!              the axis (x,y,z) coincide with the face local axes (r,s,n). In the first three columns of the array 
!              BoundaryIntegralTab() the coordinates (x,y,z) of each point are stored; in the forth column the relative d_alpha
!              (portion of solid angle relative to the point, necessary for integrations) is stored. BITcols = 4.
!              (Di Monaco et al., 2011, EACFM)                        
!----------------------------------------------------------------------------------------------------------------------------------

subroutine ComputeBoundaryIntegralTab
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: i,j,points
double precision :: delta_fi,delta_teta,fi,teta,xp,yp,zp,delta_alpha
double precision :: delta_fimz,duesendelta_fimz
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
!------------------------
! Statements
!------------------------
delta_fi = FI_INTERVAL / FI_STEPS
delta_fimz = half * delta_fi
duesendelta_fimz = two * Sin(delta_fimz)
delta_teta = TETA_INTERVAL / TETA_STEPS
points = 0
fi = -half * delta_fi
do i=1,FI_STEPS
   fi = fi + delta_fi
   teta = -half * delta_teta
   zp = -doubleh * Cos(fi)
   do j=1,TETA_STEPS
      teta = teta + delta_teta
      xp = doubleh * Sin(fi) * Cos(teta)
      yp = doubleh * Sin(fi) * Sin(teta)
      delta_alpha = delta_teta * duesendelta_fimz * Sin(fi)
      points = points + 1
      BoundIntegralTab(points,1) = xp
      BoundIntegralTab(points,2) = yp
      BoundIntegralTab(points,3) = zp
      BoundIntegralTab(points,4) = delta_alpha
! Directional cosines of the radius OP (xp,yp,zp)
      BoundIntegralTab(points,5) = xp / doubleh
      BoundIntegralTab(points,6) = yp / doubleh
      BoundIntegralTab(points,7) = zp / doubleh
   enddo
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeBoundaryIntegralTab

