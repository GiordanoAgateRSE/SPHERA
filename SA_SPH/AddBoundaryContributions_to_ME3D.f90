!cfile AddBoundaryContributions_to_ME3D.f90
!************************************************************************************
!                             S P H E R A 6.0.0
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name    : AddBoundaryContributions_to_ME3D
!
! Last updating : September 20, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07      Graphic windows calls removed
! 01  Agate/Flamini    08/10/07      Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
!AA504
! 03 Amicarelli         08/04/2014     (v5.04) Viscosity term depending on velocity divergence is removed (negligible); boundary term depending 
!                                      on molecular viscosity is disabled because of errors for granular flows and general lack of validation.
!
!************************************************************************************
! Module purpose : Module to compute boundary contributions to rodivV,
!              gradPsuro and ViscoF relative to particle npi
!              Performs implicit computation of gradPsuro
!
! Calling routine: Loop_Irre_3D
!
! Called routines: 
!
!************************************************************************************
!
  subroutine AddBoundaryContributions_to_ME3D (npi, Ncbf, tpres, tdiss, tvisc)
!
!Computes boundary contributions to rodivV, gradPsuro and ViscoF relative to particle npi
!Performs implicit computation of gradPsuro
!
!.. assign modules
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
  use FILES_ENTITIES
!
!.. Implicit Declarations ..
  implicit none
!
!.. Local Parameters ..
!  double precision,parameter :: etaratio  = 0.1d0
!  double precision,parameter :: minIntWdV = 0.001d0  
!  double precision,parameter :: pibmink = 10.0d0
!
!.. Formal Arguments ..
  integer(4),      intent(IN)    :: npi
  integer(4),      intent(IN)    :: Ncbf
  double precision,intent(INOUT),dimension(1:SPACEDIM) :: tpres, tdiss, tvisc
!
!.. Local Scalars ..
  integer(4)       :: sd, icbf, iface, ibdp
  integer(4)       :: sdj, i, j, mati, stretch
  double precision :: IntdWrm1dV
  double precision :: cinvisci, Monvisc, cinviscmult, pressi, dvn
  double precision :: FlowRate1, Lb, L, minquotanode, maxquotanode
  double precision :: Qii, roi, celeri, alfaMon, Mmult, IntGWZrm1dV
  double precision :: tpres_save1, ViscoMon_save1
!
  character(4)     :: stretchtype
!
!.. local Arrays ..
  double precision,dimension(1:SPACEDIM) :: vb, vi, dvij, RHS, nnlocal, Grav_Loc, Gpsurob_Loc, Gpsurob_Glo
  double precision,dimension(1:SPACEDIM) :: ViscoMon, ViscoShear, LocPi
!
!.. Executable Statements ..
!
!.. initializations
!
  mati = pg(npi)%imed
  cinvisci = pg(npi)%visc
  roi = pg(npi)%dens
  pressi = pg(npi)%pres
  Qii = (pressi + pressi) / roi
  tpres_save1 = zero
  ViscoMon_save1 = zero
!
  vi(:) = pg(npi)%var(:)
!
  RHS(:) = zero
  ViscoMon(:) = zero
  ViscoShear(:) = zero
!
!AA404 sub
  if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
    pg(npi)%kodvel = 0
    pg(npi)%velass(:) = zero
  end if
!
  face_loop: do icbf = 1,Ncbf
!
    ibdp = BoundaryDataPointer(3,npi) + icbf - 1 
    iface = BoundaryDataTab(ibdp)%CloBoNum
    stretch = BoundaryFace(iface)%stretch
    stretchtype = Tratto(stretch)%tipo
    nnlocal(:) = BoundaryFace(iface)%T(:,3)
!
!.. skips for the open boundary condition (no constraint must be applied in this case)
    if (stretchtype == "open") cycle face_loop
!
!!!    LocPi(:) = BoundaryDataTab(ibdp)%LocXYZ(:) !Coordinate locali della particella reale Pi
!!!    if (LocPi(3) > zero) then                  !La faccia "iface" interagisce con la particella Pi
!
    if (stretchtype == "fixe" .Or. stretchtype == "tapi") then
!
      LocPi(:) = BoundaryDataTab(ibdp)%LocXYZ(:) !Coordinate locali della particella reale Pi
      if (LocPi(3) > zero) then                  !La faccia "iface" interagisce con la particella Pi
        dvn = zero
        do SD = 1,SPACEDIM
          vb(SD) = BoundaryFace(iface)%velocity(SD)
          dvij(SD) = two * (vi(SD) - vb(SD))
          dvn = dvn + BoundaryFace(iface)%T(SD, 3) * dvij(SD)
          Grav_Loc(SD) = zero                           !Componenti locali della gravita'
          do sdj = 1,SPACEDIM
            Grav_Loc(SD) = Grav_Loc(SD) + BoundaryFace(iface)%T(sdj, SD) * Domain%grav(sdj)
          end do
        end do
           
        !*******  Boundary contribution to gradPsuro  ***********************************
            
        !Componenti locali
        do i = 1,SPACEDIM
          Gpsurob_Loc(i) = -Qii * BoundaryDataTab(ibdp)%BoundaryIntegral(3+i)
          do j = 1,SPACEDIM
            Gpsurob_Loc(i) = Gpsurob_Loc(i) - BoundaryDataTab(ibdp)%IntGiWrRdV(i, j) * Grav_Loc(j)
          end do
        end do
        !Componenti globali
        do i = 1,SPACEDIM
          Gpsurob_Glo(i) = zero
          do j = 1,SPACEDIM
            Gpsurob_Glo(i) = Gpsurob_Glo(i) + BoundaryFace(iface)%T(i, j) * Gpsurob_Loc(j)
!.. explosion
          if (esplosione) then
            tpres_save1 = tpres_save1 - (nnlocal(i) * BoundaryDataTab(ibdp)%IntGiWrRdV(i, j) &
                        + Gpsurob_Glo(i)) * dvn * nnlocal(i)
          end if
!.. explosion
          end do
          RHS(i) = RHS(i) + Gpsurob_Glo(i)
        end do
            
        !*******  Boundary contribution to ViscoF  ***********************************
        IntGWZrm1dV = BoundaryDataTab(ibdp)%BoundaryIntegral(7)
        IntdWrm1dV = BoundaryDataTab(ibdp)%BoundaryIntegral(3)
        alfaMon = Med(mati)%alfaMon
        if (alfaMon > zero .or. cinvisci > zero) then
!AA504 sub start: the molecular viscosity term, depending on velocity divergence, can be neglected (otherwise it may cause several problems); 
! further Monaghan's term is now activated even for detaching particles
            celeri = Med(mati)%celerita
            Monvisc = alfaMon * celeri * Domain%h 
!AA504 sub end
          Mmult = -Monvisc * dvn * IntGWZrm1dV
          ViscoMon(:) = ViscoMon(:) + Mmult * nnlocal(:)
!.. explosion
          do i = 1, SPACEDIM
            ViscoMon_save1 = ViscoMon_save1 + ViscoMon(i) * dvn * nnlocal(i)
          end do
!.. explosion
        end if
        if (cinvisci > zero) then
          cinviscmult = two * cinvisci * IntdWrm1dV * Tratto(stretch)%ShearCoeff
          ViscoShear(:) = ViscoShear(:) + cinviscmult * dvij(:)
        end if
      end if
!
    else if (stretchtype == "velo") then 
!
!AA404 sub
      if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%kodvel = 2
        if (tempo < Tratto(stretch)%trampa) then
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) * tempo / Tratto(stretch)%trampa
        else
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) 
        end if
      end if
      return   
! 
    else if (stretchtype == "flow") then 
!
!AA404 sub
       if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%kodvel = 2
! calcolo velocita' normale
        if (BoundaryFace(iface)%CloseParticles > 0) then
          minquotanode = 9999.0d0
          maxquotanode = const_m_9999
          do i = 1,BoundaryFace(iface)%nodes
            if (minquotanode > vertice(3,BoundaryFace(iface)%node(i)%name)) minquotanode = vertice(3,BoundaryFace(iface)%node(i)%name)
            if (maxquotanode < vertice(3,BoundaryFace(iface)%node(i)%name)) maxquotanode = vertice(3,BoundaryFace(iface)%node(i)%name)
          end do
          Lb = BoundaryFace(iface)%CloseParticles_maxQuota - minquotanode
          L = maxquotanode - minquotanode
          FlowRate1 = Tratto(stretch)%FlowRate * Lb / L
          Tratto(stretch)%NormVelocity = doubleh * FlowRate1 / (BoundaryFace(iface)%CloseParticles * Domain%PVolume)
        else
          Tratto(stretch)%NormVelocity = zero
        end if
        if (tempo < Tratto(stretch)%trampa) then
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) * tempo / Tratto(stretch)%trampa
        else
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) 
        end if
      end if
      return   
! 
!.. for the source boundary that assume three velocity components
!                                                        
    else if (stretchtype == "sour") then         
!
!AA404 sub
      if ((Domain%time_stage == 1) .or. (Domain%time_split == 1)) then 
!
        pg(npi)%kodvel = 2
        if (tempo < Tratto(stretch)%trampa) then
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) * tempo / Tratto(stretch)%trampa
        else
          pg(npi)%velass(:) = Tratto(stretch)%NormVelocity*nnlocal(:) 
        end if
      end if
      return                                                           
    end if
!
!!    end if
!
  end do face_loop

!*******  Addition of the boundary contribution to gradP  ***********************************
!.. Boundary contributions to momentum equation  
!
!!!!if (it_corrente ==1 .and. (npi== 10447 .or. npi==10425)) then
!!!!write (99,*) 'me3d',tpres,tdiss,tvisc
!!!!write (99,*) ' '
!!!!write (99,*) 'me3d',RHS,ViscoMon,ViscoShear
!!!!end if
!
  tpres(:) = tpres(:) - RHS(:)
  tdiss(:) = tdiss(:) - ViscoMon(:)
!AA504 rm and comm: this 3D boundary term has not been tested and seems not to work properly (it is erased at the moment)
! It seems useless at this stage to comment all the other lines involved as they are sparse and do not cause relevant computational time.
!  tvisc(:) = tvisc(:) - ViscoShear(:)
!............................................ 2011 mar 08
!.. contribution for Specific Internal Energy
!
  if (esplosione) then
      pg(npi)%dEdT = - half * ( tpres_save1 - ViscoMon_save1 )
  end if
!
!.. end contribution for Specific Internal Energy
!................................................
!
!prova!
!  tpres = floor(tpres * azzeramento) / azzeramento
!  where (dabs(tpres) < arrotondamento) tpres = zero
!  tdiss = floor(tdiss * azzeramento) / azzeramento
!  where (dabs(tdiss) < arrotondamento) tdiss = zero
!  tvisc = floor(tvisc * azzeramento) / azzeramento
!  where (dabs(tvisc) < arrotondamento) tvisc = zero
!prova!
!
  return
  end subroutine AddBoundaryContributions_to_ME3D
!---split

