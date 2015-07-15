!cfile ComputeBoundaryDataTab.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : ComputeBoundaryDataTab
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! 00  Agate/Guandalini  13/11/08       creation
!
!************************************************************************************
! Module purpose : Module to calculate array to store close boundary and integrals
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D
!
! Called routines: diagnostic
!                  BoundaryVolumeIntegrals2D
!                  ComputeBoundaryVolumeIntegrals_P0
!                  ComputeSurfaceIntegral_WdS2D
!                  ComputeVolumeIntegral_WdV2D
!                  FindBoundaryIntersection2D
!                  FindCloseBoundarySides2D
!                  FindCloseBoundaryFaces3D
!                  InterpolateBoundaryIntegrals2D
!                  SelectCloseBoundarySides2D
!
!************************************************************************************
!
subroutine ComputeBoundaryDataTab 
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use DIAGNOSTIC_MODULE
!
!.. Implicit Declarations ..
  implicit none
!
!.. Formal Arguments ..
!
!.. Local Scalars ..
  integer(4)       :: npi, Ncbs, IntNcbs, icbs, Ncbf, icbf, ibdt, Nfzn, Ncols
  double precision :: IntWdS, IntWdV, IntWdV1, IntdWrm1dV, IntGWZrm1dV, &
                      IntWd1s0, IntWd3s0, IntWd1s2
  double precision :: deltai, ypi, xpmin, xpmax, interlen     !xpi, 
  character(len=lencard) :: nomsub = "ComputeBoundaryDataTab"
  
!
!.. Local Arrays ..
  integer(4),      dimension(1:NUMCOLS_BIT)                   :: Colmn
  integer(4),      dimension(1:MAXCLOSEBOUNDSIDES)            :: Cloboside, Intboside
  double precision,dimension(1:PLANEDIM,1:MAXCLOSEBOUNDSIDES) :: LocXY, IntLocXY
  double precision,dimension(1:PLANEDIM)                      :: IntDpWdV
  integer(4),      dimension(1:Domain%MAXCLOSEBOUNDFACES)     :: Cloboface
  double precision,dimension(1:SPACEDIM,1:Domain%MAXCLOSEBOUNDFACES) :: LocX
  double precision,dimension(1:SPACEDIM)                      :: IntGWdV
  double precision,dimension(1:SPACEDIM,1:SPACEDIM)           :: IntGWrRdV
  double precision,dimension(1:NUMCOLS_BIT)                   :: Func
!  type (TyBoundaryData),dimension(:),allocatable :: buffer1
!
!.. Executable Statements ..
!
  BoundaryDataPointer = 0
!
  if (ncord == 2) then
!
!..  azzeramento contatore numero particelle vicine al contorno
    BoundarySide(:)%CloseParticles = 0
    BoundarySide(:)%CloseParticles_maxQuota = const_m_9999
!
!$omp parallel do default(none) &
!$omp private(npi,Ncbs,Cloboside,LocXY,IntNcbs,Intboside,IntLocXY,ibdt,icbs) &
!$omp private(xpmin,xpmax,interlen,Ncols,Colmn,deltai,Func,ypi) &
!$omp private(IntWdS,IntWdV,IntDpWdV,IntWdV1,IntWd1s0,IntWd3s0,IntWd1s2) &
!$omp shared(nag,pg,Domain,BoundaryDataTab,BoundaryDataPointer,MaxNcbs,nomsub,nout,nscr,BoundarySide,squareh)

    do npi = 1,nag
!
      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" ) cycle
!
!.. searches for the boundary sides nearest the npi-th current particle

      call FindCloseBoundarySides2D (npi,Ncbs, Cloboside, LocXY)

!
!.. some nearest boundary has been detected
      if (Ncbs > 0) then
!
!.. selects the boundary sides that effectively contribute to the motion field terms
        call SelectCloseBoundarySides2D (npi, Ncbs, Cloboside, LocXY, IntNcbs, Intboside, IntLocXY)
!
        if (IntNcbs > 0) then
!
          BoundaryDataPointer(1,npi) = Ncbs
          BoundaryDataPointer(2,npi) = IntNcbs
          ibdt = MAXCLOSEBOUNDSIDES * (npi-1)
          BoundaryDataPointer(3,npi) = ibdt+1
!
          do icbs = 1,IntNcbs
!
            ibdt = ibdt + 1
!.. controllo dimensioni vettori
            if (ibdt > MaxNcbs) then
              call diagnostic (arg1=8,arg2=1,arg3=nomsub)
!
!              MaxNcbs = ibdt * 1.2
!              write(nout,'(a,i15)') " New Max num particles*BoundaryCloseSides the new value is: MaxNcbs = ",MaxNcbs
!              allocate (buffer1(1:MaxNcbs), stat = ier)
!              if (ier /= 0) then
!                write (nout,'(1x,a,i2)') "  REDIMENSION Array BUFFER1 not allocated. Error code: ",ier
!                call diagnostic (arg1=4,arg3=nomsub)
!              else
!                write (nout,'(1x,a)') "  REDIMENSION Array BUFFER1 successfully allocated "
!              end if 
!              buffer1(1:ibdt-1) = BoundaryDataTab(1:ibdt-1)
!              deallocate (BoundaryDataTab)
!              allocate (BoundaryDataTab(1:MaxNcbs), stat = ier)
!              if (ier /= 0) then
!                write (nout,'(1x,a,i2)') "  REDIMENSION Array BoundaryDataTab not allocated. Error code: ",ier
!                call diagnostic (arg1=4,arg3=nomsub)
!              else
!                write (nout,'(1x,a)') "  REDIMENSION Array BoundaryDataTab successfully allocated "
!              end if 
!              BoundaryDataTab(1:ibdt-1) = buffer1(1:ibdt-1)
!              deallocate (buffer1)
            end if
!
            BoundaryDataTab(ibdt)%CloBoNum  = Intboside(icbs)
            BoundaryDataTab(ibdt)%LocXYZ(1:PLANEDIM) = IntLocXY(1:PLANEDIM,icbs)
            BoundaryDataTab(ibdt)%LocXYZ(3) = zero
            Func = zero
!
            if (IntNcbs == 2) then  !Calcolo numerico degli integrali IntWds e IntWdV
!.. Calcolo integrali
!.. 2D computation of the boundary integrals (semy-analitic approach, volumetric method)
!.. to find the intersection between the boundary and the kernel support
              call FindBoundaryIntersection2D (icbs, Intboside, IntLocXY, BoundarySide, xpmin, xpmax, interlen)
!.. computation of the 1D integrals
              call ComputeSurfaceIntegral_WdS2D (icbs, IntLocXY, xpmin, interlen, IntWdS)  !Intboside, xpmax, 
!.. computation of the 2D integrals
              call BoundaryVolumeIntegrals2D (icbs, IntLocXY, xpmin, xpmax, interlen, IntWdV, IntDpWdV)
              call ComputeVolumeIntegral_WdV2D (icbs, IntNcbs, Intboside, IntLocXY, BoundarySide, xpmin, xpmax, interlen, IntWdV1)
!.. Interpolation of some integrals (from tables)
!.. Interpolazione da tabella degli integrali IntWd1s0, IntWd3s0, IntWd1s2 
              Ncols=3
              Colmn(1)=3
              Colmn(2)=4
              Colmn(3)=5
              ypi = IntLocXY(2, icbs)
              deltai = ypi / Domain%h
              call InterpolateBoundaryIntegrals2D (deltai, Ncols, Colmn, Func)
              IntWd1s0 = Func(1) / squareh
              IntWd3s0 = Func(2) / squareh
              IntWd1s2 = Func(3) / squareh
!
              BoundaryDataTab(ibdt)%BoundaryIntegral(1) = IntWdS
              BoundaryDataTab(ibdt)%BoundaryIntegral(2) = IntWdV
              BoundaryDataTab(ibdt)%BoundaryIntegral(3) = IntWdV1
              BoundaryDataTab(ibdt)%BoundaryIntegral(4:5) = IntDpWdV(1:2)
              BoundaryDataTab(ibdt)%BoundaryIntegral(6) = IntWd1s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(7) = IntWd3s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(8) = IntWd1s2
! 
            else if (IntNcbs == 1) then
 ! Interpolazione da tabella degli integrali IntWdS,IntWdV, IntWd1s0, IntWd3s0, IntWd1s2 
              Ncols=5
              Colmn(1) = 1
              Colmn(2) = 2
              Colmn(3) = 3
              Colmn(4) = 4
              Colmn(5) = 5
!              xpi = IntLocXY(1, icbs)
              ypi = IntLocXY(2, icbs)
              deltai = ypi / Domain%h
              call InterpolateBoundaryIntegrals2D (deltai, Ncols, Colmn, Func)
              IntWdS = Func(1) / Domain%h
              IntWdV = Func(2)
              IntWd1s0 = Func(3) / squareh
              IntWd3s0 = Func(4) / squareh
              IntWd1s2 = Func(5) / squareh
              BoundaryDataTab(ibdt)%BoundaryIntegral(1) = IntWds
              BoundaryDataTab(ibdt)%BoundaryIntegral(2) = IntWdV
              BoundaryDataTab(ibdt)%BoundaryIntegral(3) = IntWdV
              BoundaryDataTab(ibdt)%BoundaryIntegral(4:5) = zero
              BoundaryDataTab(ibdt)%BoundaryIntegral(6) = IntWd1s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(7) = IntWd3s0
              BoundaryDataTab(ibdt)%BoundaryIntegral(8) = IntWd1s2
            end if
          end do
        end if
      end if
    end do
!
!$omp end parallel do
!
!----------------------
!.. 3D
!----------------------
  else
!
!..  azzeramento contatore numero particelle vicine al contorno
    BoundaryFace(:)%CloseParticles = 0
    BoundaryFace(:)%CloseParticles_maxQuota = const_m_9999
!
    
!AA504sub
!$omp parallel do default(none) &
!$omp private(npi,Ncbf,Cloboface,LocX,Nfzn,icbf,ibdt) &
!$omp private(IntWdV,IntdWrm1dV,IntGWZrm1dV,IntGWdV,IntGWrRdV) &
!$omp shared(nag,pg,BoundaryDataTab,BoundaryDataPointer,EpCount,MaxNcbf,nout,nscr,nomsub,Domain)

    loop_particle:  do npi = 1,nag
!
      if (pg(npi)%cella == 0 .or. pg(npi)%vel_type /= "std" ) cycle loop_particle
!
!.. searches for the boundary faces nearest the npi-th current particle

      call FindCloseBoundaryFaces3D (npi, Ncbf, Cloboface, LocX, Nfzn)

!
      if (Ncbf == 0) then
!
        BoundaryDataPointer(:,npi) = 0
!
      else
!
!.. some nearest boundary has been detected
!
!.. conteggio particelle sfuggite dal contorno ma ancora all'interno della griglia (EpCount)
!..  coordinate normali alle faces tutte negative
        if (Nfzn == Ncbf) then
          EpCount(pg(npi)%imed) = EpCount(pg(npi)%imed) + 1
!.. eliminazione della particella
!          write (699,'(a,f10.4,a,i10,a,3f15.6)') 'time ',tempo,' escaped particle n.',npi,' coordinates ', pg(npi)%coord
!          pg(nag)%cella = -2
        end if
!
        BoundaryDataPointer(1,npi) = Ncbf
        BoundaryDataPointer(2,npi) = 0
        ibdt = Domain%MAXCLOSEBOUNDFACES * (npi-1)
        BoundaryDataPointer(3,npi) = ibdt + 1
!
!.. controllo dimensioni vettori
        if (ibdt > MaxNcbf) then
          call diagnostic (arg1=8,arg2=2,arg3=nomsub)
        end if
!
        do icbf = 1,Ncbf
!
          call ComputeBoundaryVolumeIntegrals_P0 (icbf, Cloboface, LocX, IntWdV, IntdWrm1dV, IntGWZrm1dV, IntGWdV, IntGWrRdV)
!
          ibdt = ibdt + 1
!
          BoundaryDataTab(ibdt)%CloBoNum  = CloboFace(icbf)
          BoundaryDataTab(ibdt)%LocXYZ(1:SPACEDIM)    = LocX(1:SPACEDIM,icbf)
          BoundaryDataTab(ibdt)%BoundaryIntegral(1)   = zero
          BoundaryDataTab(ibdt)%BoundaryIntegral(2)   = IntWdV
          BoundaryDataTab(ibdt)%BoundaryIntegral(3)   = IntdWrm1dV
          BoundaryDataTab(ibdt)%BoundaryIntegral(4:6) = IntGWdV(1:3)
          BoundaryDataTab(ibdt)%BoundaryIntegral(7)   = IntGWZrm1dV
          BoundaryDataTab(ibdt)%BoundaryIntegral(8)   = zero
          BoundaryDataTab(ibdt)%IntGiWrRdV(:,:)       = IntGWrRdV(:,:)
!
        end do
!
      end if
!
    end do loop_particle
!
!AA504 sub
!$omp end parallel do
!
  end if
!
return
end subroutine ComputeBoundaryDataTab
!---split

