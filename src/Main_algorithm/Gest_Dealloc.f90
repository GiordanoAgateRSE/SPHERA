!----------------------------------------------------------------------------------------------------------------------------------
! SPHERA (Smoothed Particle Hydrodynamics research software; mesh-less Computational Fluid Dynamics code).
! Copyright 2005-2015 (RSE SpA -formerly ERSE SpA, formerly CESI RICERCA, formerly CESI-) 
!      
!     
!   
!      
!  

! This file is part of SPHERA.
!  
!  
!  
!  
! SPHERA is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
!  
!  
!  
!----------------------------------------------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------------------------------------------
! Program unit: Gest_Dealloc         
! Description: Deallocations.                   
!----------------------------------------------------------------------------------------------------------------------------------

subroutine Gest_Dealloc(nomsub)
!------------------------
! Modules
!------------------------ 
use I_O_file_module
use Static_allocation_module
use Dynamic_allocation_module
use I_O_diagnostic_module
!------------------------
! Declarations
!------------------------
implicit none
logical :: check
integer(4) :: ier,i
character(LEN=lencard), intent(IN) :: nomsub
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
write (nout,'(1x,a)') ">> Storage deallocation in routine "//trim(nomsub)
check = .true.
if (allocated(Vertice)) then
   deallocate(Vertice,stat=ier)
      if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: VERTICE not deallocated with error code: ",ier
         check = .false.
      else
         write (nout,'(1x,a)') "   Array: VERTICE successfully deallocated "
   endif
endif
if (allocated(Tratto)) then
   deallocate(Tratto,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: TRATTO not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: TRATTO successfully deallocated "
   endif
endif
if (allocated(BoundaryFace)) then
   deallocate(BoundaryFace,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: BoundaryFace not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: BoundaryFace successfully deallocated"
   endif
endif
if (allocated(BFaceList)) then 
   deallocate(BFaceList,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: BFACELIST not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: BFACELIST successfully deallocated "
   endif
endif 
if (allocated(BoundaryVertex)) then
   deallocate(BoundaryVertex,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: BOUNDARYVERTEX not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: BOUNDARYVERTEX successfully dealloc."
   endif
endif
if (allocated(BoundarySide)) then
   deallocate(BoundarySide,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: BOUNDARYSIDE not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: BOUNDARYSIDE successfully deallocated"
   endif
endif
if (allocated(Partz)) then
   deallocate(Partz,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: PARTZ not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: PARTZ successfully deallocated "
   endif
endif
if (allocated(Med)) then
   deallocate(Med,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
      "   Array: MED not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: MED successfully deallocated "
   endif
endif
if (allocated(OpCount)) then
   deallocate(OpCount,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: OpCount not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: OpCount successfully deallocated "
   endif
endif
if (allocated(SpCount)) then
   deallocate(SpCount,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: SpCount not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: SpCount successfully deallocated "
   endif
endif
if (allocated(EpCount)) then
   deallocate(EpCount,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: EpCount not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: EpCount successfully deallocated "
   endif
endif
if (allocated(EpOrdGrid)) then
   deallocate(EpOrdGrid,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: EpOrdGrid not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: EpOrdGrid successfully deallocated "
   endif
endif
if (allocated(Control_Sections)) then
   deallocate(Control_Sections,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: CONTROL_SECTION not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: CONTROL_SECTION successfully deallocated"
   endif
endif
if (allocated(Control_Points)) then
   deallocate(Control_Points,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: CONTROL_POINTS not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: CONTROL_POINTS successfully deallocated"
   endif
endif
if (allocated(Control_Lines)) then
   deallocate(Control_Lines,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: CONTROL_LINES not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: CONTROL_LINES successfully deallocated "
   endif
endif
if (allocated(Pg)) then
   deallocate(Pg,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: PG not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: PG successfully deallocated "
   endif
endif
if (allocated(ts0_pg)) then
   deallocate(ts0_pg,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
      "   Array: ts0_pg not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: ts0_pg successfully deallocated "
   endif
endif
if (allocated(Section_Points)) then
   deallocate(Section_Points,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: SECTION_POINTS not deallocated with error code: ",ier
      check = .false.
      else
      write (nout,'(1x,a)') "   Array: SECTION_POINTS successfully deallocated "
   end if
endif
if (allocated(NPartOrd)) then
   deallocate(NPartOrd,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: NPARTORD not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: NPARTORD successfully deallocated "
   endif
endif
if (allocated(Icont)) then
   deallocate(Icont,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: ICONT not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: ICONT successfully deallocated "
   endif
endif
if (allocated(nPartIntorno)) then
   deallocate(nPartIntorno,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: nPartIntorno not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: nPartIntorno successfully deallocated"
   endif
endif
if (allocated(PartIntorno)) then
   deallocate(PartIntorno,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: PartIntorno not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: PartIntorno successfully deallocated "
   endif
endif
if (allocated(PartKernel)) then
   deallocate(PartKernel,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: PartKernel not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: PartKernel successfully deallocated "
   endif
endif
if (allocated(rag)) then
   deallocate(rag,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: rag not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: rag successfully deallocated "
   endif
endif
if (allocated(BoundaryDataTab)) then
   deallocate(BoundaryDataTab,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 & 
         "   Array: BoundaryDataTab not deallocated with error code: ",ier
      check = .false.
      else
      write (nout,'(1x,a)') "   Array: BoundaryDataTab successfully deallocated"
   endif
endif
if (allocated(BoundaryDataPointer)) then
   deallocate(BoundaryDataPointer,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: BoundaryDataPointer not deallocated with error code: ",ier
      check = .false.
      else
      write (nout,'(1x,a)')                                                    &
         "   Array: BoundaryDataPointer successfully deallocated"
   endif
endif
if (allocated(Array_Flu)) then
   deallocate(Array_Flu,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: Array_Flu not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: Array_Flu successfully deallocated "
   endif
endif
if (Domain%tipo=="bsph") then
   if (allocated(NPartOrd_w)) then
      deallocate(NPartOrd_w,stat=ier)
      if (ier/=0) then
         write (nout,'(1x,a,i2)')                                              &
            "   Array: NPARTORD_w not deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)')                                              &
               "   Array: NPARTORD_w successfully deallocated "
      endif
   endif
   if (allocated(Icont_w)) then
      deallocate(Icont_w,stat=ier)
      if (ier/=0) then
         write (nout,'(1x,a,i2)')                                              &
            "   Array: ICONT_w not deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)') "   Array: ICONT_w successfully deallocated "
      endif
   endif
   if (allocated(nPartIntorno_fw)) then
      deallocate(nPartIntorno_fw,stat=ier)
      if (ier/=0) then
         write (nout,'(1x,a,i2)')                                              &
         "   Array: nPartIntorno_fw not deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)')                                              &
               "   Array: nPartIntorno_fw successfully deallocated "
      endif
   endif
   if (allocated(PartIntorno_fw)) then
      deallocate(PartIntorno_fw,stat=ier)
      if (ier/=0) then
         write (nout,'(1x,a,i2)')                                              &
            "   Array PartIntorno_fw not deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)')                                              &
               "   Array PartIntorno_fw successfully deallocated "
      endif
   endif
   if (allocated(grad_vel_VSL_fw)) then
      deallocate(grad_vel_VSL_fw,stat=ier)
      if (ier/=0) then
         write (nout,'(1x,a,i2)')                                              &
            "   Array grad_vel_VSL_fw not deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)')                                              &
               "   Array grad_vel_VSL_fw successfully deallocated "
      endif
   endif
   if (allocated(kernel_fw)) then
      deallocate(kernel_fw,stat=ier)
      if (ier/=0) then
         write (nout,'(1x,a,i2)')                                              &
            "   Array kernel_fw not deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)') "   Array kernel_fw successfully deallocated"
        endif
   endif
   if (allocated(rag_fw)) then
      deallocate(rag_fw,stat=ier)
      if (ier/=0) then
         write (nout,'(1x,a,i2)')                                              &
            "   Array: rag_fw not deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)') "   Array: rag_fw successfully deallocated "
      endif
   endif
endif
do i=1,n_bodies
   if (allocated(body_arr(i)%body_kinematics)) then
      deallocate(body_arr(i)%body_kinematics,stat =ier)
      if (ier/=0) then
         write (nout,'(1x,2a,i2)') "   Array: body_arr(i)%body_kinematics not",&
             " deallocated with error code: ",ier
         check = .false.
         else
            write (nout,'(1x,a)')                                              &
               "   Array: body_arr(i)%body_kinematics successfully deallocated "
      endif
   endif
enddo
if (allocated(body_arr)) then
   deallocate(body_arr,stat =ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: body_arr not deallocated with error code: ",ier
      check = .false.
      else
      write (nout,'(1x,a)') "   Array: body_arr successfully deallocated "
   endif
endif
if (allocated(bp_arr)) then
   deallocate(bp_arr,stat =ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: bp_arr not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: bp_arr successfully deallocated "
   endif
endif
if (allocated(Icont_bp)) then
   deallocate(Icont_bp,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: Icont_bp not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: Icont_bp successfully deallocated "
   endif
endif
if (allocated(NPartOrd_bp)) then
   deallocate(NPartOrd_bp,stat =ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: NPartOrd_bp not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: NPartOrd_bp successfully deallocated "
   endif
endif
if (allocated(nPartIntorno_bp_f)) then
   deallocate(nPartIntorno_bp_f,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
      "   Array: nPartIntorno_bp_f not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: nPartIntorno_bp_f successfully deallocated "
   endif
endif
if (allocated(PartIntorno_bp_f)) then
   deallocate(PartIntorno_bp_f,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: PartIntorno_bp_f not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: PartIntorno_bp_f successfully deallocated "
   endif
endif
if (allocated(KerDer_bp_f_cub_spl)) then
   deallocate(KerDer_bp_f_cub_spl,stat =ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: KerDer_bp_f_cub_spl not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: KerDer_bp_f_cub_spl successfully deallocated "
   endif
endif
if (allocated(KerDer_bp_f_Gal)) then
   deallocate(KerDer_bp_f_Gal,stat =ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: KerDer_bp_f_Gal not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: KerDer_bp_f_Gal successfully deallocated "
   endif
endif 
if (allocated(rag_bp_f)) then
   deallocate(rag_bp_f,stat =ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: rag_bp_f not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: rag_bp_f successfully deallocated "
   endif
endif
if (allocated(surf_body_part)) then
   deallocate(surf_body_part,stat =ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: surf_body_part not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: surf_body_part successfully deallocated "
   endif
endif 
if (allocated(nPartIntorno_bp_bp)) then
   deallocate(nPartIntorno_bp_bp,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: nPartIntorno_bp_bp not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: nPartIntorno_bp_bp successfully deallocated "
   endif
endif
if (allocated(PartIntorno_bp_bp)) then
   deallocate(PartIntorno_bp_bp,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: PartIntorno_bp_bp not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)')                                                 &
            "   Array: PartIntorno_bp_bp successfully deallocated "
   endif
endif
if (allocated(rag_bp_bp)) then
   deallocate(rag_bp_bp,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: rag_bp_bp not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: rag_bp_bp successfully deallocated "
   endif
endif
if (allocated(impact_vel)) then
   deallocate(impact_vel,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: impact_vel not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: impact_vel successfully deallocated "
   endif
endif
if (allocated(BoundaryConvexEdge)) then
   deallocate(BoundaryConvexEdge,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: BoundaryConvexEdge not deallocated with error code: ",ier
      check = .false.
      else
      write (nout,'(1x,a)')                                                    &
         "   Array: BoundaryConvexEdge successfully deallocated "
   endif
endif
if (allocated(GCBFVector)) then
   deallocate(GCBFVector,stat=ier)
   if (ier/=0) then
      write (nout,'(1x,a,i2)')                                                 &
         "   Array: GCBFVector not deallocated with error code: ",ier
      check = .false.
      else
         write (nout,'(1x,a)') "   Array: GCBFVector successfully deallocated "
   endif
endif 
if(allocated(Q_sections%section)) then
   do i=1,Q_sections%n_sect
      if (allocated(Q_sections%section(i)%flow_rate))                          &
         deallocate(Q_sections%section(i)%flow_rate) 
   enddo    
   deallocate(Q_sections%section)
endif 
if(allocated(Granular_flows_options%lines)) then
   deallocate(Granular_flows_options%lines)
endif
if (allocated(DBSPH%kinematics)) then
   deallocate(DBSPH%kinematics,STAT=ier)
   if (ier/=0) then
      write(nout,*) 'Deallocation of aux_array in GestDealloc failed; ',       &
         'the program terminates here.'
      call diagnostic (arg1=5,arg2=340)       
      else
         write (nout,'(1x,a)')                                                 &
            "Deallocation of aux_array in GestDealloc successfully completed."
   endif
endif
if (allocated(DBSPH%inlet_sections)) then
   deallocate(DBSPH%inlet_sections,STAT=ier)
   if (ier/=0) then
      write(nout,*) 'Deallocation of DBSPH%inlet_sections in GestDealloc ',    &
         'failed; the program terminates here.'
         call diagnostic (arg1=5,arg2=340)       
      else
         write (nout,'(1x,2a)') "Deallocation of DBSPH%inlet_sections in ",    &
            "GestDealloc successfully completed."
   endif
endif
if (allocated(DBSPH%outlet_sections)) then
   deallocate(DBSPH%outlet_sections,STAT=ier)
   if (ier/=0) then
      write(nout,*) 'Deallocation of DBSPH%outlet_sections in GestDealloc ',   &
         'failed; the program terminates here.'
      call diagnostic (arg1=5,arg2=340)       
      else
         write (nout,'(1x,2a)') "Deallocation of DBSPH%outlet_sections in ",   &
            "GestDealloc successfully completed."
   endif
endif
if (check) return
call diagnostic(arg1=3,arg3=nomsub)
!------------------------
! Deallocations
!------------------------
return
end subroutine Gest_dealloc

