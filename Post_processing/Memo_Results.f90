!cfile Memo_Results.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Memo_Results
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!
!************************************************************************************
! Module purpose : Module to write results for post-processor
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, diagnostic
!
! Called routines: 
!
!************************************************************************************
!
subroutine Memo_Results ( it, it_memo, it_rest, dtvel, str)
!
!.. assign modules
use FILES_ENTITIES
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
implicit none
!
!.. Formal Arguments ..
integer(4),      intent(IN)    :: it
integer(4),      intent(INOUT) :: it_memo
integer(4),      intent(INOUT) :: it_rest
double precision,intent(IN)    :: dtvel
character(6),    intent(IN)    :: str
!
!.. Local Scalars ..
integer(4) :: nrecords, restartcode
!
!
!.. Executable Statements ..
!
restartcode = 0
!
!=== SCRITTURA SU FILE RISULTATI e di RESTART ######################
!
! caso step a scrittura forzata (iniziale o finale)
! if ( it == it_memo ) then
 if ( index(str,'inizio')/=0 ) then
!
! memorizzazione informazioni fisse
    nrecords = 5
    if ( NumVertici   > 0 ) nrecords = nrecords + 1
    if ( NumFacce     > 0 ) nrecords = nrecords + 1
    if ( NumFacce     > 0 ) nrecords = nrecords + 1
    if ( NumTratti    > 0 ) nrecords = nrecords + 1
    if ( NPartZone    > 0 ) nrecords = nrecords + 1
    if ( NumBVertices > 0 ) nrecords = nrecords + 1
    if ( NumBSides    > 0 ) nrecords = nrecords + 1
!    if ( NPointst     > 0 ) nrecords = nrecords + 1
!    if ( NLines       > 0 ) nrecords = nrecords + 1
!    if ( NSections    > 0 ) nrecords = nrecords + 1
    write (nres) version,nrecords
    write (nres) Ncord, Nag, NMedium, NPartZone, &
                 NumVertici, NumFacce, NumTratti, NumBVertices, NumBSides, &
                 NPointst,NPoints,NPointsl,NPointse, NLines, NSections, GCBFVecDim, &
                 doubleh
    write (nres) domain
    write (nres) grid
    write (nres) Med(1:NMedium)
    if ( NumVertici   > 0 ) write (nres) Vertice(1:SPACEDIM,1:NumVertici)
    if ( NumFacce     > 0 ) write (nres) BoundaryFace(1:NumFacce)
    if ( NumFacce     > 0 ) write (nres) BFaceList(1:NumFacce)
    if ( NumTratti    > 0 ) write (nres) Tratto(1:NumTratti)
    if ( NPartZone    > 0 ) write (nres) Partz(1:NPartZone)
    if ( NumBVertices > 0 ) write (nres) BoundaryVertex(1:NumBVertices)
    if ( NumBSides    > 0 ) write (nres) BoundarySide(1:NumBSides)
!    if ( NPointst     > 0 ) write (nres) control_points(1:NPointst)
!    if ( NLines       > 0 ) write (nres) control_lines(1:NLines)
!    if ( NSections    > 0 ) write (nres) Control_Sections(0:NSections+1)
    flush(nres)
!
    write (nout,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
    write (nout,'(a,i10,a,f15.5)') " Results and restart heading saved   step: ",it,"   time: ",tempo
    write (nout,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
    write (nscr,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
    write (nscr,'(a,i10,a,f15.5)') " Results and restart heading saved   step: ",it,"   time: ",tempo
    write (nscr,'(a,i10,a,f15.5)') " ----------------------------------------------------------------------------"
!
  end if
!  else 
! memorizzazione informazioni variabili nel tempo
! caso memorizzazione nel loop
! 
! caso delta step con restart 
    if (Domain%irest_fr > 0) then
     if ( mod(it,Domain%irest_fr) == 0 ) then
       it_rest = it
     end if
! caso delta time con restart 
    else if (Domain%rest_fr > zero) then
     if ( it > 1 .and. mod(tempo,Domain%rest_fr) <= dtvel ) then
       it_rest = it
     end if
    end if
!
! caso delta step con solo risultati 
    if (Domain%imemo_fr > 0) then
     if ( mod(it,Domain%imemo_fr) == 0 ) then
        it_memo = it
      end if
! caso delta time con solo risultati 
     else if (Domain%memo_fr > zero) then
      if ( it > 1 .and. mod(tempo,Domain%memo_fr) <= dtvel ) then
        it_memo = it
      end if
    end if
!
!    if (it_rest == it) then
    if (it_rest == it .or. index(str,'inizio') /= 0 .or. index(str,'fine') /= 0) then
! memorizzo per restart
      restartcode = 1 ! se restartcode=1 memorizzo tutto il vettore pg
      write (nres) it,tempo,dt,nag,ncord,restartcode
      write (nres) pg(1:nag)
      flush(nres)
!
      if ( index(str,'inizio') == 0 ) then
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
        write (nout,'(a,i10,a,f15.5)') " Results and restart saved   step: ",it,"   time: ",tempo
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " Results and restart saved   step: ",it,"   time: ",tempo
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------------------"
      end if
!
    else if (it_memo == it) then
! memorizzo per risultati
      restartcode = 0 ! se restartcode=0 memorizzo il vettore pg solo per visualizzazione 
      write (nres) it,tempo,dt,nag,ncord,restartcode
      write (nres) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3), &
                   pg(1:nag)%vel(1),  pg(1:nag)%vel(2),  pg(1:nag)%vel(3), &
                   pg(1:nag)%pres,   &
                   pg(1:nag)%dens,   &
                   pg(1:nag)%mass,   &
                   pg(1:nag)%visc,   &
                   pg(1:nag)%IntEn,  &
                   pg(1:nag)%VolFra, &
                   pg(1:nag)%imed,   &
                   pg(1:nag)%icol
      flush(nres)
!
      if ( index(str,'inizio') == 0 ) then
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------"
        write (nout,'(a,i10,a,f15.5)') " Results saved   step: ",it,"   time: ",tempo
        write (nout,'(a,i10,a,f15.5)') " --------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------"
        write (nscr,'(a,i10,a,f15.5)') " Results saved   step: ",it,"   time: ",tempo
        write (nscr,'(a,i10,a,f15.5)') " --------------------------------------------------------"
      end if
    end if
!
! end if
!
  return
  end subroutine Memo_Results
!---split

