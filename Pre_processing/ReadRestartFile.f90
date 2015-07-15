!cfile ReadRestartFile.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : readrestartfile
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
! Module purpose : Module to read results file for rstart purpose
!
! Calling routine: gest_input
!
! Called routines: diagnostic
!
!************************************************************************************

subroutine ReadRestartFile ( option, ier, nrecords )
!
!.. assign modules
use files_entities
use GLOBAL_MODULE
use AdM_USER_TYPE
use ALLOC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
character(7),intent(IN)  :: option
integer(4),intent(INOUT) :: ier
integer(4),intent(INOUT) :: nrecords 
!
!.. Local Scalars ..
!
integer(4)       :: restartcode, save_istart
integer(4)       :: ioerr
double precision :: save_start
character(12)    :: ainp = "Restart File"
character(len=8) :: versionerest
!
!.. External Routines ..
character(80), external :: lcase
logical,       external :: ReadCheck
!
!.. Executable Statements ..
!
 ier = 0
!
!.. Restart heading 
!
 if ( TRIM(lcase(option)) == TRIM(lcase("heading")) ) then
!
   rewind (nsav)
!
   write(nout,'(a)')    "-------------------"
   write(nout,"(1x,a)") ">> Restart heading."
   write(nout,'(a)')    "-------------------"
   write(nscr,'(a)')    "-------------------"
   write(nscr,"(1x,a)") ">> Restart heading."
   write(nscr,'(a)')    "-------------------"
!
   read (nsav,iostat=ioerr) versionerest,nrecords
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"versionerest,nrecords",nsav,nout) ) return
! check sulla versione
   if (TRIM(lcase(version)) /= TRIM(lcase(versionerest))) then
     write(nout,'(a)')    "---------------------------------------------------------------"
     write(nout,"(1x,a)") ">> ERROR! The Restart version is not equal the current version."
     write(nout,"(1x,a)") ">>        The Run is stopped."
     write(nout,'(a)')    "---------------------------------------------------------------"
     flush(nout)
     write(nscr,'(a)')    "---------------------------------------------------------------"
     write(nscr,"(1x,a)") ">> ERROR! The Restart version is not equal the current version."
     write(nscr,"(1x,a)") ">>        The Run is stopped."
     write(nscr,'(a)')    "---------------------------------------------------------------"
     flush(nscr)
     stop
   end if
!
   read (nsav,iostat=ioerr) Ncord, Nag, NMedium, NPartZone, &
                            NumVertici, NumFacce, NumTratti, NumBVertices, NumBSides, &
                            NPointst,NPoints,NPointsl,NPointse, NLines, NSections, GCBFVecDim, &
                            doubleh
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"ncord, nag, ...",nsav,nout) ) return
!
 else if ( TRIM(lcase(option)) == "reading" ) then
!
   write(nout,'(a)')    "-----------------------------------------------------------------------"
   write(nout,"(1x,a)") ">> Restart reading:  step          time      interval    num.particles"
   write(nscr,'(a)')    "-----------------------------------------------------------------------"
   write(nscr,"(1x,a)") ">> Restart reading:  step          time      interval    num.particles"
!
   save_istart = Domain%istart
   save_start  = Domain%start
!
   read (nsav,iostat=ioerr) domain
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"domain",nsav,nout) ) return
!
   read (nsav,iostat=ioerr) grid
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"grid",nsav,nout) ) return
!.. allocazione matrice 2d per calcolo pelolibero (caso erosione)
!AA504 sub (removed everywhere the fifth element of the array
  allocate (ind_interfaces(Grid%ncd(1),Grid%ncd(2),4), stat = ioerr)
  if (ioerr /= 0) then
!AA504 sub      
    write (nout,'(1x,a,i2)') "    Array ind_interfaces not allocated. Error code: ",ioerr
    stop ' routine ReadRestartFile'
  else
!AA504 sub      
    write (nout,'(1x,a)') "    Array ind_interfaces successfully allocated "
  end if
!
   read (nsav,iostat=ioerr) Med(1:NMedium)
   if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Med",nsav,nout) ) return
!
   if ( NumVertici   > 0 ) then
     read (nsav,iostat=ioerr) Vertice(1:SPACEDIM,1:NumVertici)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Vertice",nsav,nout) ) return
   end if
!
   if ( NumFacce     > 0 ) then 
     read (nsav,iostat=ioerr) BoundaryFace(1:NumFacce)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BoundaryFace",nsav,nout) ) return
   end if
!
   if ( NumFacce     > 0 ) then
     read (nsav,iostat=ioerr) BFaceList(1:NumFacce)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BFaceList",nsav,nout) ) return
   end if
!
   if ( NumTratti    > 0 ) then
     read (nsav,iostat=ioerr) Tratto(1:NumTratti)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Tratto",nsav,nout) ) return
   end if
!
   if ( NPartZone    > 0 ) then
     read (nsav,iostat=ioerr) Partz(1:NPartZone)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"Partz",nsav,nout) ) return
   end if
!
   if ( NumBVertices > 0 ) then
     read (nsav,iostat=ioerr) BoundaryVertex(1:NumBVertices)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BoundaryVertex",nsav,nout) ) return
   end if
!
   if ( NumBSides    > 0 ) then
     read (nsav,iostat=ioerr) BoundarySide(1:NumBSides)
     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"BoundarySide",nsav,nout) ) return
   end if
!
!   if ( NPointst     > 0 ) then
!     read (nsav,iostat=ioerr) control_points(1:NPointst)
!     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"control_points",nsav,nout) ) return
!   end if
!
!   if ( NLines       > 0 ) then
!     read (nsav,iostat=ioerr) control_lines(1:NLines)
!     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"control_lines",nsav,nout) ) return
!   end if
!
!   if ( NSections    > 0 ) then
!     read (nsav,iostat=ioerr) Control_Sections(0:NSections+1)
!     if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"control_sections",nsav,nout) ) return
!   end if
!
!
!.. Restart positioning is based on the step number
!
   it_start = 0 
!
!   if ( Domain%istart > 0 ) then
   if ( save_istart > 0 ) then
!
!     do while ( Domain%istart > it_start )
     do while ( save_istart > it_start )
!
       read(nsav,iostat=ioerr) it_start,tempo,dt,nag,ncord,restartcode
       if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"it_start,tempo,dt,nag,ncord,restartcode",nsav,nout) ) return
!
       write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nout)
       write(nscr,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nscr)
!
!       if ( it_start < Domain%istart ) then
       if ( it_start < save_istart ) then
         read (nsav,iostat=ioerr) 
         if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"...",nsav,nout) ) return
!
       else
! leggo per restart
         if (restartcode == 1) then
           read (nsav,iostat=ioerr) pg(1:nag)
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
!
           write(nout,'(a)') " "
           write(nout,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nout)
           write(nscr,'(a)') " "
           write(nscr,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nscr)
! leggo per risultati
         else if (restartcode == 0) then
           read (nsav,iostat=ioerr) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3), &
                                    pg(1:nag)%vel(1),  pg(1:nag)%vel(2),  pg(1:nag)%vel(3), &
                                    pg(1:nag)%pres,   &
                                    pg(1:nag)%dens,   &
                                    pg(1:nag)%mass,   &
                                    pg(1:nag)%visc,   &
                                    pg(1:nag)%IntEn,  &
                                    pg(1:nag)%VolFra, &
                                    pg(1:nag)%imed,   &
                                    pg(1:nag)%icol
!
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
!
           write(nout,'(a)') " "
           write(nout,'(a,i10,a,g12.5)') "   Located Result Step :",it_start,"   Time :",tempo; flush(nout)
           write(nout,'(a)') "       But this step is not a restart step. Check the correct step for restart in the restart file."; flush(nout)
           write(nout,'(a)') " The program is terminated."; flush(nout)
           write(nscr,'(a)') " "
           write(nscr,'(a,i10,a,g12.5)') "   Located Result Step :",it_start,"   Time :",tempo; flush(nscr)
           write(nscr,'(a)') "       But this step is not a restart step. Check the correct step for restart in the restart file."; flush(nscr)
           write(nscr,'(a)') " The program is terminated."; flush(nscr)
           stop
         end if
         return
       end if
     end do
     write(nout,'(a,i10,a)') "   Restart Step Number:",it_start," has not been found"
     write(nscr,'(a,i10,a)') "   Restart Step Number:",it_start," has not been found"
     ier = 3
!
!.. Restart positioning is based on the time step
!
   else if ( save_start > zero ) then
!
     tempo = zero
!
     do while ( save_start > tempo )
!
       read(nsav,iostat=ioerr) it_start,tempo,dt,nag,ncord,restartcode
       if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"it_start,tempo,dt,nag,ncord,restartcode",nsav,nout) ) return
!
       write(nout,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nout)
       write(nscr,"(16x,i10,2(2x,g12.5),7x,i10)") it_start,tempo,dt,nag; flush(nscr)
!
       if ( tempo < Domain%start ) then
         read (nsav,iostat=ioerr) 
         if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"...",nsav,nout) ) return
       else
! leggo per restart
         if (restartcode == 1) then
           read (nsav,iostat=ioerr) pg(1:nag)
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
           write(nout,'(a)') 
           write(nout,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nout)
           write(nscr,'(a)') 
           write(nscr,'(a,i10,a,g12.5)') "   Located Restart Step :",it_start,"   Time :",tempo; flush(nscr)
! leggo per risultati
         else if (restartcode == 0) then
           read (nsav,iostat=ioerr) pg(1:nag)%coord(1),pg(1:nag)%coord(2),pg(1:nag)%coord(3), &
                                    pg(1:nag)%vel(1),  pg(1:nag)%vel(2),  pg(1:nag)%vel  (3), &
                                    pg(1:nag)%pres,   &
                                    pg(1:nag)%dens,   &
                                    pg(1:nag)%mass,   &
                                    pg(1:nag)%visc,   &
                                    pg(1:nag)%IntEn,  &
                                    pg(1:nag)%VolFra, &
                                    pg(1:nag)%imed,   &
                                    pg(1:nag)%icol
           if ( .NOT.ReadCheck (ioerr,ier,it_start,ainp,"pg",nsav,nout) ) return
           write(nout,'(a)') 
           write(nout,'(a,i10,a,g12.5)') "   Located Result Time :",it_start,"   Time :",tempo; flush(nout)
           write(nout,'(a)') "       But this time is not a restart time. Check the correct time for restart in the restart file."; flush(nout)
           write(nout,'(a)') " The program is terminated."; flush(nout)
           write(nscr,'(a)') 
           write(nscr,'(a,i10,a,g12.5)') "   Located Result Time :",it_start,"   Time :",tempo; flush(nscr)
           write(nscr,'(a)') "       But this time is not a restart time. Check the correct time for restart in the restart file."; flush(nscr)
           write(nscr,'(a)') " The program is terminated."; flush(nscr)
           stop
         end if
         return
       end if
     end do
     write(nout,'(a,i10,a)') "   Restart Time Step:",Domain%start," has not been found"
     write(nscr,'(a,i10,a)') "   Restart Time Step:",Domain%start," has not been found"
     ier = 3
!
   else
     write (nout,'(a)' ) "  > Restart cannot be read at step:",it_start,"  time:",tempo
     write (nscr,'(a)' ) "  > Restart cannot be read at step:",it_start,"  time:",tempo
     ier = 4
   end if
!
   write (nout,'(a)' ) "  > Restart read successfully at step:",it_start,"  time:",tempo
   write (nscr,'(a)' ) "  > Restart read successfully at step:",it_start,"  time:",tempo
!
 else
!
   ier = 5
!
 end if
!
return
endsubroutine ReadRestartFile
!---split

