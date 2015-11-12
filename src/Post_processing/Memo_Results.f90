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
! Program unit: Memo_Results               
! Description: To write detailed results for restart. Not recommended.       
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Memo_Results(it,it_memo,it_rest,dtvel,str)
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Hybrid_allocation_module
use Dynamic_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4),intent(IN) :: it
double precision,intent(IN) :: dtvel
character(6),intent(IN) :: str
integer(4),intent(INOUT) :: it_memo
integer(4),intent(INOUT) :: it_rest
integer(4) :: nrecords, restartcode
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
restartcode = 0
!------------------------
! Statements
!------------------------
! Post-processing for first (and last) step 
if (index(str,'inizio')/=0) then
   nrecords = 5
   if (NumVertici>0) nrecords = nrecords + 1
   if (NumFacce>0) nrecords = nrecords + 1
   if (NumFacce>0) nrecords = nrecords + 1
   if (NumTratti>0) nrecords = nrecords + 1
   if (NPartZone>0) nrecords = nrecords + 1
   if (NumBVertices>0) nrecords = nrecords + 1
   if (NumBSides>0) nrecords = nrecords + 1
   write(nres) nrecords
   write(nres) Ncord,Nag,NMedium,NPartZone,NumVertici, NumFacce, NumTratti,    &
      NumBVertices,NumBSides,NPointst,NPoints,NPointsl,NPointse,NLines,        &
      NSections, GCBFVecDim,doubleh
   write(nres) domain
   write(nres) grid
   write(nres) Med(1:NMedium)
   if (NumVertici>0) write(nres) Vertice(1:SPACEDIM,1:NumVertici)
   if (NumFacce>0) write(nres) BoundaryFace(1:NumFacce)
   if (NumFacce>0) write(nres) BFaceList(1:NumFacce)
   if (NumTratti>0) write(nres) Tratto(1:NumTratti)
   if (NPartZone>0) write(nres) Partz(1:NPartZone)
   if (NumBVertices>0) write(nres) BoundaryVertex(1:NumBVertices)
   if (NumBSides>0) write(nres) BoundarySide(1:NumBSides)
   flush(nres)
   write(nout,'(a,i10,a,f15.5)')                                               &
" ----------------------------------------------------------------------------"
   write(nout,'(a,i10,a,f15.5)') " Results and restart heading saved   step: ",&
      it,"   time: ",tempo
   write(nout,'(a,i10,a,f15.5)')                                               &
" ----------------------------------------------------------------------------"
   write(nscr,'(a,i10,a,f15.5)')                                               &
" ----------------------------------------------------------------------------"
   write(nscr,'(a,i10,a,f15.5)') " Results and restart heading saved   step: ",&
      it,"   time: ",tempo
   write(nscr,'(a,i10,a,f15.5)')                                               &
" ----------------------------------------------------------------------------"
endif
if (Domain%irest_fr>0) then
   if (mod(it,Domain%irest_fr)==0) then
      it_rest = it
   endif
! Case with restart  
   elseif (Domain%rest_fr>zero) then
      if ((it>1).and.(mod(tempo,Domain%rest_fr)<=dtvel)) then
         it_rest = it
     endif
endif
if (Domain%imemo_fr>0) then
   if (mod(it,Domain%imemo_fr)==0) then
        it_memo = it
   endif
   elseif (Domain%memo_fr>zero) then
      if ((it>1).and.(mod(tempo,Domain%memo_fr)<=dtvel)) then
         it_memo = it
      endif
endif
if ((it_rest==it).or.(index(str,'inizio')/=0).or.(index(str,'fine')/=0)) then
! If restartcode=1, then to save the whole array "pg"
   restartcode = 1 
   write(nres) it,tempo,dt,nag,ncord,restartcode
   write(nres) pg(1:nag)
   flush(nres)
   if (index(str,'inizio')==0) then
      write(nout,'(a,i10,a,f15.5)')                                            &
      " --------------------------------------------------------------------"
      write(nout,'(a,i10,a,f15.5)') " Results and restart saved   step: ",it,  &
         "   time: ",tempo
      write(nout,'(a,i10,a,f15.5)')                                            &
      " --------------------------------------------------------------------"
      write(nscr,'(a,i10,a,f15.5)')                                            &
      " --------------------------------------------------------------------"
      write(nscr,'(a,i10,a,f15.5)') " Results and restart saved   step: ",it,  &
         "   time: ",tempo
      write(nscr,'(a,i10,a,f15.5)')                                            &
      " --------------------------------------------------------------------"
   endif
   elseif (it_memo==it) then
! If restartcode=0, then to save "pg" only for visualizations
      restartcode = 0  
      write(nres) it,tempo,dt,nag,ncord,restartcode
      write(nres) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3),    &
         pg(1:nag)%vel(1),pg(1:nag)%vel(2),pg(1:nag)%vel(3),pg(1:nag)%pres,    &
         pg(1:nag)%dens,pg(1:nag)%mass,pg(1:nag)%visc,pg(1:nag)%IntEn,         &
         pg(1:nag)%VolFra,pg(1:nag)%imed,pg(1:nag)%icol
      flush(nres)
      if (index(str,'inizio')==0) then
         write(nout,'(a,i10,a,f15.5)')                                         &
            " --------------------------------------------------------"
         write(nout,'(a,i10,a,f15.5)') " Results saved   step: ",it,"   time: "&
            ,tempo
         write(nout,'(a,i10,a,f15.5)')                                         &
            " --------------------------------------------------------"
         write(nscr,'(a,i10,a,f15.5)')                                         &
            " --------------------------------------------------------"
         write(nscr,'(a,i10,a,f15.5)') " Results saved   step: ",it,"   time: "&
            ,tempo
         write(nscr,'(a,i10,a,f15.5)')                                         &
            " --------------------------------------------------------"
      endif
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine Memo_Results

