!-------------------------------------------------------------------------------
! SPHERA v.8.0 (Smoothed Particle Hydrodynamics research software; mesh-less
! Computational Fluid Dynamics code).
! Copyright 2005-2017 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA,
! formerly CESI-Ricerca di Sistema)
!
! SPHERA authors and email contact are provided on SPHERA documentation.
!
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
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Program unit: ComputeBoundaryVolumeIntegrals_P0                                   
! Description: (Di Monaco et al., 2011, EACFM)                        
!-------------------------------------------------------------------------------
subroutine ComputeBoundaryVolumeIntegrals_P0(icbf,Clobface,LocX,IntWdV,        &
   IntdWrm1dV,IntGWZrm1dV,IntGWdV,IntGWrRdV)
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
integer(4),intent(IN) :: icbf
integer(4),intent(IN),dimension(1:Domain%MAXCLOSEBOUNDFACES) :: Clobface
double precision,intent(IN),dimension(1:SPACEDIM,1:Domain%MAXCLOSEBOUNDFACES) :: LocX
double precision,intent(OUT) :: IntWdV,IntdWrm1dV,IntGWZrm1dV
double precision,intent(OUT),dimension(1:SPACEDIM) :: IntGWdV
double precision,intent(OUT),dimension(1:SPACEDIM,1:SPACEDIM) :: IntGWrRdV
integer(4),parameter :: nicols = 4
double precision,parameter :: eps = 0.05d0
integer(4) :: SD,sdj,ipk,fkod,iface
double precision :: deltaPiPx,PiPxdist,rob,delta_alpha,tsegnato,zpmin,RZ,RZ2 
double precision :: robpre,JW3ro2dA,JdW3ro1dA,JdW3ro2dA,JdW3ro3dA
integer(4),dimension(1:nicols) :: icol
double precision,dimension(1:SPACEDIM) :: LocPi,LocPj,LocPiPj,PXLoc,csi,RLoccos
double precision,dimension(1:nicols) :: ivalue
logical, external :: IsPointInternal
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
icol(1) = 1
icol(2) = 2
icol(3) = 3
icol(4) = 4
IntWdV = zero
IntdWrm1dV = zero
IntGWZrm1dV = zero
IntGWdV(:) = zero
IntGWrRdV(:,:) = zero
JW3ro2dA  = zero
JdW3ro1dA = zero
JdW3ro2dA = zero
JdW3ro3dA = zero
iface = Clobface(icbf)
zpmin = eps * Domain%h
!------------------------
! Statements
!------------------------
! Local coordinates of particle Pi
LocPi(:) = LocX(:, icbf)                     
if (LocPi(3)<zero) return
if (LocPi(3)<zpmin) LocPi(3) = zpmin
robpre = 2.1
do ipk=1,BITrows
! Local components of the vector "PiPj"
   LocPiPj(:) = BoundIntegralTab(ipk,1:3)      
! Local coordinates of point Pj
   LocPj(:) = LocPi(:) + LocPiPj(:)             
   if (LocPj(3)<=0) then   
! Point Pj does not lie on the same semi-space of Pi, 
! with respect to the face "iface"
! Local coorinates of point "PX", which is the intersection of the segment 
! "PiPj" with the face "iface" 
      tsegnato = LocPi(3) / (LocPi(3) - LocPj(3))
      PXLoc(:) = LocPi(:) + (LocPj(:) - LocPi(:)) * tsegnato
      PXLoc(3) = zero
      call LocalNormalCoordinates (PXLoc, csi, iface)
      fkod = BoundaryFace(iface)%nodes - 2
      if (IsPointInternal(fkod, csi)) then    
! "PX" belongs to the face. Thus, Pj contributes to the boundary integral
! Distance between "Pi" and "Px" 
         PiPxdist = zero
         do SD=1,SPACEDIM
            deltaPiPx = PXLoc(SD) - LocPi(SD)
            PiPxdist = PiPxdist + deltaPiPx * deltaPiPx
         enddo
         PiPxdist = Dsqrt(PiPxdist)
! Normalised distance and related functions 
         delta_alpha = BoundIntegralTab(ipk, 4)
         rob = PiPxdist / Domain%h
         if (Abs(rob - robpre)>0.001) then
            ivalue = zero
            call InterpolateTable(rob,nicols,icol,ivalue)
            JW3ro2dA  = ivalue(1) * delta_alpha
            JdW3ro1dA = ivalue(2) * delta_alpha
            JdW3ro2dA = ivalue(3) * delta_alpha
            JdW3ro3dA = ivalue(4) * delta_alpha
            robpre = rob
         endif
         RLoccos(:) = BoundIntegralTab(ipk, 5:7)
         RZ = RLoccos(3)
         RZ2 = RZ * RZ
! Boundary volume integrals  
         IntWdV = IntWdV + JW3ro2dA
         IntdWrm1dV = IntdWrm1dV + JdW3ro1dA
         IntGWZrm1dV = IntGWZrm1dV + RZ2 * JdW3ro1dA
         do SD=1,SPACEDIM
            IntGWdV(SD) = IntGWdV(SD) + RLoccos(SD) * JdW3ro2dA
            do sdj=SD,SPACEDIM
               IntGWrRdV(SD,sdj) = IntGWrRdV(SD,sdj) + RLoccos(SD) *           &
                  RLoccos(sdj) * JdW3ro3dA
            enddo
         enddo
      endif
   endif
enddo
! Completing the symmetric matrix "IntGWrRdV(3,3)"
IntGWrRdV(2,1) = IntGWrRdV(1,2)
IntGWrRdV(3,1) = IntGWrRdV(1,3)
IntGWrRdV(3,2) = IntGWrRdV(2,3)
if (LocPi(3)==zpmin) then
   PXLoc(:) = LocPi(:)
   PXLoc(3) = zero
   call LocalNormalCoordinates (PXLoc, csi, iface)
   fkod = BoundaryFace(iface)%nodes - 2
   if (IsPointInternal(fkod, csi)) then    
! The particle projection belongs to the face
      IntWdV = half
      else                                    
! The particle projection is external to the face
         IntWdV = zero
   endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine ComputeBoundaryVolumeIntegrals_P0

