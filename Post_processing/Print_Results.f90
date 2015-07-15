!cfile Print_Results.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : Print_Results
!
! Last updating : November 23, 2011
!
! Improvement traceback:
!
! ..  E.Bon, A. Di Monaco, S. Falappi  Initial development of the code
! 00  Agate/Guandalini  28/08/07       Graphic windows calls removed
! 01  Agate/Flamini     08/10/07       Check of entire code
! 02  Agate/Guandalini  2008           Check and review entire code
! 03  Agate             2011           varie
! 05  Amicarelli        23/11/2011     multiple inlet
! 06  Amicarelli/Agate  30Nov11        BSPH: wall element parameters
!AA501b comment
! 07  Amicarelli-Agate  13nov12        Body dynamics
!AA504
! 08  Amicarelli        08apr14        (v5.04) Modifications for 3D erosion criterion and elapsed time
!
!************************************************************************************
! Module purpose : Module to write results on output ascii file
!
! Calling routine: Loop_Irre_2D, Loop_Irre_3D, diagnostic
!
! Called routines: 
!
!************************************************************************************
!
  subroutine Print_Results (it, it_print, str)
!
!.. assign modules
  use FILES_ENTITIES
  use GLOBAL_MODULE
  use AdM_USER_TYPE
  use ALLOC_MODULE
!
!.. implicit declarations
  implicit none
!
!.. dummy arguments
  integer(4),  intent(IN)    :: it
  integer(4),  intent(INOUT) :: it_print
  character(6),intent(IN)    :: str
!
!.. local Scalars ..
  integer(4)        :: npi,i,codice,dummy,OpCountot,SpCountot
  integer(4)        :: minlocvelo,maxlocvelo
  integer(4)        :: minlocvelx,maxlocvelx
  integer(4)        :: minlocvely,maxlocvely
  integer(4)        :: minlocvelz,maxlocvelz
  integer(4)        :: minlocpres,maxlocpres
  integer(4)        :: minlocdens,maxlocdens
  integer(4)        :: minlocvisc,maxlocvisc
  integer(4)        :: minloccodi,maxloccodi
  integer(4)        :: minlocInEn,maxlocInEn
  double precision  :: minvelo,maxvelo
  double precision  :: minvelx,maxvelx
  double precision  :: minvely,maxvely
  double precision  :: minvelz,maxvelz
  double precision  :: minpres,maxpres
  double precision  :: mindens,maxdens
  double precision  :: minvisc,maxvisc
  double precision  :: mincodi,maxcodi
  double precision  :: minInEn,maxInEn
  double precision  :: modvel
  character(len=42) :: fmt100="(a,i10,a,e18.9,a,e18.9,a,i  ,a,i  ,a,i  ))"
  character(len=47) :: fmt101="(a,2(1x,f10.4,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=47) :: fmt102="(a,2(1x,f10.1,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=47) :: fmt103="(a,2(1x,f10.6,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=47) :: fmt104="(a,2(1x,g10.4,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
!AA504  
  character(len=47) :: fmt105="(a,2(1x,e10.4,1x,a,1x,i6,1x,a,3(1x,f8.2,1x,a)))"
  character(len=12) :: stringa
  character(len=2)  :: coppia
!
!AA406 start
  integer(4)        :: minlocvelo_w,maxlocvelo_w,minlocpres_w,maxlocpres_w
  double precision  :: minvelo_w,maxvelo_w,minpres_w,maxpres_w
!AA406 end

!AA501b start
  integer(4)        :: minlocvelo_bp,maxlocvelo_bp,minlocpres_bp,maxlocpres_bp
  double precision  :: minvelo_bp,maxvelo_bp,minpres_bp,maxpres_bp
  integer(4)        :: minlocvelo_body,maxlocvelo_body,minlocomega_body,maxlocomega_body,nbi
  double precision  :: minvelo_body,maxvelo_body,minomega_body,maxomega_body,modomega
!AA501b end

!AA504 start
  integer(4)        :: minloctau_tauc,maxloctau_tauc,minlock_BetaGamma,maxlock_BetaGamma,minlocu_star,maxlocu_star
  integer(4)        :: machine_Julian_day,machine_hour,machine_minute,machine_second
  double precision  :: mintau_tauc,maxtau_tauc,mink_BetaGamma,maxk_BetaGamma,minu_star,maxu_star,time_elapsed_tot_est
!AA504 end

!
!.. local Arrays ..
  integer(4),  dimension(1) :: pos
!
!.. executable statements
!
!.. detects the output code for prints
!
  codice = abs(Domain%ioutopt)
  if ( index(str,'inizio')/=0 .or. index(str,'fine')/=0) then
    if ( codice == 0 ) return
  else
    if ( codice == 0 .or. mod(it,codice) /= 0 ) return
  end if
!
!.. evaluates the proper formats
!
  dummy = count(pg(1:nag)%cella >0)
  write(stringa,'(i12)') dummy
  i = len_trim(adjustl(stringa))
  write (coppia,'(i2)') i
  fmt100(27:28) = adjustl(coppia)
!
  OpCountot = 0
  SpCountot = 0
  do i = 1,NMedium
    OpCountot = OpCountot + OpCount(i)
    SpCountot = SpCountot + SpCount(i)
  end do
  write(stringa,'(i12)') OpCountot
  i = len_trim(adjustl(stringa))
  write (coppia,'(i2)') i
  fmt100(33:34) = adjustl(coppia)
!
  write(stringa,'(i12)') SpCountot
  i = len_trim(adjustl(stringa))
  write (coppia,'(i2)') i
  fmt100(39:40) = adjustl(coppia)
!
!AA406 
  if (nag>0) then
!
!.. searches for the minimum and maximum values of all the distributions and 
!.. the corresponding particle identifiers
!
!.. module of the velocity and their components
!
  minvelo = max_positive_number
  maxvelo = max_negative_number
!
!AA406 start
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
     minvelo_w = max_positive_number
     maxvelo_w = max_negative_number
  endif
!AA406 end
!

!AA501b start
 if (n_bodies > 0) then
     minvelo_bp = max_positive_number
     maxvelo_bp = max_negative_number
     minvelo_body = max_positive_number
     maxvelo_body = max_negative_number
     minomega_body = max_positive_number
     maxomega_body = max_negative_number
  endif
!AA501b end

!AA504 start
 if ((Granular_flows_options%erosion_flag.ne.1).and.(Granular_flows_options%ID_erosion_criterion==1)) then
    mintau_tauc = max_positive_number
    maxtau_tauc = max_negative_number
    minu_star = max_positive_number
    maxu_star = max_negative_number
    if (Granular_flows_options%ID_erosion_criterion == 1) then    
       mink_BetaGamma = max_positive_number
       maxk_BetaGamma = max_negative_number
    endif
 endif
!AA504 end

  do npi = 1,nag
!   
    if (pg(npi)%cella == 0) cycle
    modvel = pg(npi)%vel(1)*pg(npi)%vel(1)+pg(npi)%vel(2)*pg(npi)%vel(2)+pg(npi)%vel(3)*pg(npi)%vel(3)
    if (modvel < minvelo) then
      minvelo = modvel
      minlocvelo = npi
    end if
    if (modvel > maxvelo) then
      maxvelo = modvel
      maxlocvelo = npi
    end if
 !
  end do
  minvelo = Dsqrt(minvelo)
  maxvelo = Dsqrt(maxvelo)
!
  minvelx = minval(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  maxvelx = maxval(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  minlocvelx = pos(1)
  pos = maxloc(pg(1:nag)%vel(1),mask=pg(1:nag)%cella /= 0)
  maxlocvelx = pos(1)
!
  minvely = minval(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  maxvely = maxval(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  minlocvely = pos(1)
  pos = maxloc(pg(1:nag)%vel(2),mask=pg(1:nag)%cella /= 0)
  maxlocvely = pos(1)
!
  minvelz = minval(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  maxvelz = maxval(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  minlocvelz = pos(1)
  pos = maxloc(pg(1:nag)%vel(3),mask=pg(1:nag)%cella /= 0)
  maxlocvelz = pos(1)
!
!.. pressure
!
  minpres = minval(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  maxpres = maxval(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  minlocpres = pos(1)
  pos = maxloc(pg(1:nag)%pres,mask=pg(1:nag)%cella /= 0)
  maxlocpres = pos(1)
!
!.. density
!
  mindens = minval(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  maxdens = maxval(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  minlocdens = pos(1)
  pos = maxloc(pg(1:nag)%dens,mask=pg(1:nag)%cella /= 0)
  maxlocdens = pos(1)
!
!.. viscosity
!
  minvisc = minval(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  maxvisc = maxval(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  pos = minloc(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  minlocvisc = pos(1)
  pos = maxloc(pg(1:nag)%visc,mask=pg(1:nag)%cella /= 0)
  maxlocvisc = pos(1)
!
!.. diffusion coefficient
!
  if (diffusione) then
    mincodi = minval(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    maxcodi = maxval(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    pos = minloc(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    minloccodi = pos(1)
    pos = maxloc(pg(1:nag)%coefdif,mask=pg(1:nag)%cella /= 0)
    maxloccodi = pos(1)
  end if
!
!.. esplosion coefficient
!
  if (esplosione) then
    minInEn = minval(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    maxInEn = maxval(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    pos = minloc(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    minlocInEn = pos(1)
    pos = maxloc(pg(1:nag)%IntEn,mask=pg(1:nag)%cella /= 0)
    maxlocInEn = pos(1)
  end if
!
!
!AA406 start
!Wall parameter limits
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
! Wall pressure 
     minpres_w = minval(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     maxpres_w = maxval(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     pos = minloc(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     minlocpres_w = pos(1)
     pos = maxloc(pg_w(1:DBSPH%n_w)%pres,mask=pg_w(1:DBSPH%n_w)%cella /= 0)
     maxlocpres_w = pos(1)
! Wall velocity
     do npi = 1,DBSPH%n_w
        if (pg_w(npi)%cella == 0) cycle
        modvel = pg_w(npi)%vel(1)*pg_w(npi)%vel(1)+pg_w(npi)%vel(2)*pg_w(npi)%vel(2)+pg_w(npi)%vel(3)*pg_w(npi)%vel(3)
        if (modvel < minvelo_w) then
           minvelo_w = modvel
           minlocvelo_w = npi
        end if
        if (modvel > maxvelo_w) then
           maxvelo_w = modvel
           maxlocvelo_w = npi
        end if
     end do
     minvelo_w = Dsqrt(minvelo_w)
     maxvelo_w = Dsqrt(maxvelo_w)
  endif
!AA406 end

!AA501b start
!limits for body particles and bodies
 if (n_bodies > 0) then
! pressure (body particles)
     minpres_bp = minval(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     maxpres_bp = maxval(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     pos = minloc(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     minlocpres_bp = pos(1)
     pos = maxloc(bp_arr(1:n_body_part)%pres,mask=bp_arr(1:n_body_part)%cell /= 0)
     maxlocpres_bp = pos(1)
! velocity (body particles)
     do npi = 1,n_body_part
        if (bp_arr(npi)%cell == 0) cycle
        modvel = bp_arr(npi)%vel(1)*bp_arr(npi)%vel(1)+bp_arr(npi)%vel(2)*bp_arr(npi)%vel(2)+bp_arr(npi)%vel(3)*bp_arr(npi)%vel(3)
        if (modvel < minvelo_bp) then
           minvelo_bp = modvel
           minlocvelo_bp = npi
        end if
        if (modvel > maxvelo_bp) then
           maxvelo_bp = modvel
           maxlocvelo_bp = npi
        end if
     end do
     minvelo_bp = Dsqrt(minvelo_bp)
     maxvelo_bp = Dsqrt(maxvelo_bp)
! velocity (bodies)
     do nbi = 1,n_bodies
        modvel = body_arr(nbi)%u_CM(1)*body_arr(nbi)%u_CM(1)+body_arr(nbi)%u_CM(2)*body_arr(nbi)%u_CM(2)+ &
                 body_arr(nbi)%u_CM(3)*body_arr(nbi)%u_CM(3)
        if (modvel < minvelo_body) then
           minvelo_body = modvel
           minlocvelo_body = nbi
        end if
        if (modvel > maxvelo_body) then
           maxvelo_body = modvel
           maxlocvelo_body = nbi
        end if
     end do
     minvelo_body = Dsqrt(minvelo_body)
     maxvelo_body = Dsqrt(maxvelo_body)
! angular velocity (bodies)
     do nbi = 1,n_bodies
        modomega = dsqrt(dot_product(body_arr(nbi)%omega,body_arr(nbi)%omega))
        if (modomega < minomega_body) then
           minomega_body = modomega
           minlocomega_body = nbi
        end if
        if (modomega > maxomega_body) then
           maxomega_body = modomega
           maxlocomega_body = nbi
        end if
     end do
     minomega_body = Dsqrt(minomega_body)
     maxomega_body = Dsqrt(maxomega_body)
  endif
!AA501b end

!AA504 start
!Limits for supplementary bed load transport parameters 
 if ((Granular_flows_options%erosion_flag.ne.1).and.(Granular_flows_options%ID_erosion_criterion==1)) then
! tau_tauc 
     mintau_tauc = minval(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     maxtau_tauc = maxval(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     pos = minloc(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     minloctau_tauc = pos(1)
     pos = maxloc(pg(1:nag)%tau_tauc,mask=pg(1:nag)%cella /= 0)
     maxloctau_tauc = pos(1)
! u_star  
     minu_star = minval(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     maxu_star = maxval(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     pos = minloc(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     minlocu_star = pos(1)
     pos = maxloc(pg(1:nag)%u_star,mask=pg(1:nag)%cella /= 0)
     maxlocu_star = pos(1) 
! k_BetaGamma 
     if (Granular_flows_options%ID_erosion_criterion == 1) then   
        mink_BetaGamma = minval(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        maxk_BetaGamma = maxval(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        pos = minloc(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        minlock_BetaGamma = pos(1)
        pos = maxloc(pg(1:nag)%k_BetaGamma,mask=pg(1:nag)%cella /= 0)
        maxlock_BetaGamma = pos(1)      
     endif   
 endif
!AA504 end

!AA504 start
 if (exetype == "linux") then
    if (Domain%tmax>0.d0) then
       call system("date +%j%H%M%S > date_now.txt")
       open (unit_time_elapsed,file='date_now.txt',status="unknown",form="formatted")
       read (unit_time_elapsed,'(i3,i2,i2,i2)') machine_Julian_day,machine_hour,machine_minute,machine_second
       close (unit_time_elapsed)
       call system("rm -f date_now.txt")
    endif   
 endif 
!AA504 end
 
!
!.. final prints
!
!AA406 rm (moved some lines before)
!AA405 
!  if (nag>0) then
!
  write (nout,'(128("."))') 
  write (nout,fmt100) &
  " Print at:     | step: ",it," | time: ",tempo," | Dt: ",dt," | Particles: inside ",dummy," gone out ",OpCountot, &
  " gone in ",SpCountot
!AA504 start  
  if (exetype == "linux") then
    time_elapsed_tot_est = ((Domain%t_pre_iter-Domain%t0) + &
                            (machine_Julian_day*24*60*60+machine_hour*60*60+machine_minute*60+machine_second - Domain%t_pre_iter) * 1.d0) / (3600.0d0)  
    if (time_elapsed_tot_est<0.d0) time_elapsed_tot_est = time_elapsed_tot_est + 366.d0*24.d0*60.d0*60.d0   
    write (nout,'(a,g12.5,a,g12.5,a)') "Elapsed time: ",time_elapsed_tot_est," hours = ",time_elapsed_tot_est/24.d0," days."
    time_elapsed_tot_est = ((Domain%t_pre_iter-Domain%t0) + &
                            (machine_Julian_day*24*60*60+machine_hour*60*60+machine_minute*60+machine_second - Domain%t_pre_iter) &
                             * (Domain%tmax/tempo)) / (3600.0d0)  
    if (time_elapsed_tot_est<0.d0) time_elapsed_tot_est = time_elapsed_tot_est + 366.d0*24.d0*60.d0*60.d0  
    write (nout,'(a,g12.5,a,g12.5,a)') "Elapsed time (at the end of the simulation, real time estimation): ",time_elapsed_tot_est," hours = ", &
                                        time_elapsed_tot_est/24.d0," days."
  end if
!AA504 end  
  write (nout,fmt101) &
  " ............. |  Min. val. |Min.loc.| X coord. | Y coord. | Z coord. ||  Max. val. |Max.loc.| X coord. | Y coord. | Z coord. |"
  write (nout,fmt101) &
  "  Tot velocity |",minvelo,"|",minlocvelo,"|",pg(minlocvelo)%coord(1),"|",pg(minlocvelo)%coord(2),"|",pg(minlocvelo)%coord(3), &
  "||",maxvelo,"|",maxlocvelo,"|",pg(maxlocvelo)%coord(1),"|",pg(maxlocvelo)%coord(2),"|",pg(maxlocvelo)%coord(3),"|"
  write (nout,fmt101) &
  "  velocity x   |",minvelx,"|",minlocvelx,"|",pg(minlocvelx)%coord(1),"|",pg(minlocvelx)%coord(2),"|",pg(minlocvelx)%coord(3), &
  "||",maxvelx,"|",maxlocvelx,"|",pg(maxlocvelx)%coord(1),"|",pg(maxlocvelx)%coord(2),"|",pg(maxlocvelx)%coord(3),"|"
  write (nout,fmt101) &
  "  velocity y   |",minvely,"|",minlocvely,"|",pg(minlocvely)%coord(1),"|",pg(minlocvely)%coord(2),"|",pg(minlocvely)%coord(3), &
  "||",maxvely,"|",maxlocvely,"|",pg(maxlocvely)%coord(1),"|",pg(maxlocvely)%coord(2),"|",pg(maxlocvely)%coord(3),"|"
  write (nout,fmt101) &
  "  velocity z   |",minvelz,"|",minlocvelz,"|",pg(minlocvelz)%coord(1),"|",pg(minlocvelz)%coord(2),"|",pg(minlocvelz)%coord(3), &
  "||",maxvelz,"|",maxlocvelz,"|",pg(maxlocvelz)%coord(1),"|",pg(maxlocvelz)%coord(2),"|",pg(maxlocvelz)%coord(3),"|"
  if (esplosione) then
    write (nout,fmt104) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nout,fmt104) &
  "  Int.Energy   |",minInEn,"|",minlocInEn,"|",pg(minlocInEn)%coord(1),"|",pg(minlocInEn)%coord(2),"|", &
  pg(minlocInEn)%coord(3),"||",maxInEn,"|",maxlocInEn,"|",pg(maxlocInEn)%coord(1),"|",pg(maxlocInEn)%coord(2),"|", &
  pg(maxlocInEn)%coord(3),"|"
    write (nout,fmt104) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  else
    write (nout,fmt102) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nout,fmt102) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  end if
!AA504 sub  
  write (nout,fmt105) &
  "  viscosity    |",minvisc,"|",minlocvisc,"|",pg(minlocvisc)%coord(1),"|",pg(minlocvisc)%coord(2),"|",pg(minlocvisc)%coord(3), &
  "||",maxvisc,"|",maxlocvisc,"|",pg(maxlocvisc)%coord(1),"|",pg(maxlocvisc)%coord(2),"|",pg(maxlocvisc)%coord(3),"|"
  if (diffusione) then
    write (nout,fmt103) &
  "  coef.diff.   |",mincodi,"|",minloccodi,"|",pg(minloccodi)%coord(1),"|",pg(minloccodi)%coord(2),"|", &
  pg(minloccodi)%coord(3),"||",maxcodi,"|",maxloccodi,"|",pg(maxloccodi)%coord(1),"|",pg(maxloccodi)%coord(2),"|", &
  pg(maxloccodi)%coord(3),"|"
  end if
  write (nscr,'(128("."))') 
  write (nscr,fmt100) &
  " Print at:     | step: ",it," | time: ",tempo," | Dt: ",dt," | Particles: inside ",dummy," gone out ",OpCountot, &
  " gone in ",SpCountot
  write (nscr,fmt101) &
  " ............. |  Min. val. |Min.loc.| X coord. | Y coord. | Z coord. ||  Max. val. |Max.loc.| X coord. | Y coord. | Z coord. |"
  write (nscr,fmt101) &
  "  Tot velocity |",minvelo,"|",minlocvelo,"|",pg(minlocvelo)%coord(1),"|",pg(minlocvelo)%coord(2),"|",pg(minlocvelo)%coord(3), &
  "||",maxvelo,"|",maxlocvelo,"|",pg(maxlocvelo)%coord(1),"|",pg(maxlocvelo)%coord(2),"|",pg(maxlocvelo)%coord(3),"|"
  write (nscr,fmt101) &
  "  velocity x   |",minvelx,"|",minlocvelx,"|",pg(minlocvelx)%coord(1),"|",pg(minlocvelx)%coord(2),"|",pg(minlocvelx)%coord(3), &
  "||",maxvelx,"|",maxlocvelx,"|",pg(maxlocvelx)%coord(1),"|",pg(maxlocvelx)%coord(2),"|",pg(maxlocvelx)%coord(3),"|"
  write (nscr,fmt101) &
  "  velocity y   |",minvely,"|",minlocvely,"|",pg(minlocvely)%coord(1),"|",pg(minlocvely)%coord(2),"|",pg(minlocvely)%coord(3), &
  "||",maxvely,"|",maxlocvely,"|",pg(maxlocvely)%coord(1),"|",pg(maxlocvely)%coord(2),"|",pg(maxlocvely)%coord(3),"|"
  write (nscr,fmt101) &
  "  velocity z   |",minvelz,"|",minlocvelz,"|",pg(minlocvelz)%coord(1),"|",pg(minlocvelz)%coord(2),"|",pg(minlocvelz)%coord(3), &
  "||",maxvelz,"|",maxlocvelz,"|",pg(maxlocvelz)%coord(1),"|",pg(maxlocvelz)%coord(2),"|",pg(maxlocvelz)%coord(3),"|"
  if (esplosione) then
    write (nscr,fmt104) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nscr,fmt104) &
  "  Int.Energy   |",minInEn,"|",minlocInEn,"|",pg(minlocInEn)%coord(1),"|",pg(minlocInEn)%coord(2),"|",&
  pg(minlocInEn)%coord(3),"||",maxInEn,"|",maxlocInEn,"|",pg(maxlocInEn)%coord(1),"|",pg(maxlocInEn)%coord(2),"|", &
  pg(maxlocInEn)%coord(3),"|"
    write (nscr,fmt104) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  else
    write (nscr,fmt102) &
  "  pressure     |",minpres,"|",minlocpres,"|",pg(minlocpres)%coord(1),"|",pg(minlocpres)%coord(2),"|",pg(minlocpres)%coord(3), &
  "||",maxpres,"|",maxlocpres,"|",pg(maxlocpres)%coord(1),"|",pg(maxlocpres)%coord(2),"|",pg(maxlocpres)%coord(3),"|"
    write (nscr,fmt101) &
  "  density      |",mindens,"|",minlocdens,"|",pg(minlocdens)%coord(1),"|",pg(minlocdens)%coord(2),"|",pg(minlocdens)%coord(3), &
  "||",maxdens,"|",maxlocdens,"|",pg(maxlocdens)%coord(1),"|",pg(maxlocdens)%coord(2),"|",pg(maxlocdens)%coord(3),"|"
  end if
!AA504sub
  write (nscr,fmt105) &
  "  viscosity    |",minvisc,"|",minlocvisc,"|",pg(minlocvisc)%coord(1),"|",pg(minlocvisc)%coord(2),"|",pg(minlocvisc)%coord(3), &
  "||",maxvisc,"|",maxlocvisc,"|",pg(maxlocvisc)%coord(1),"|",pg(maxlocvisc)%coord(2),"|",pg(maxlocvisc)%coord(3),"|"
  if (diffusione) then
    write (nscr,fmt103) &
  "  coef.diff.   |",mincodi,"|",minloccodi,"|",pg(minloccodi)%coord(1),"|",pg(minloccodi)%coord(2),"|",&
  pg(minloccodi)%coord(3),"||",maxcodi,"|",maxloccodi,"|",pg(maxloccodi)%coord(1),"|",pg(maxloccodi)%coord(2),"|", &
  pg(maxloccodi)%coord(3),"|"
  end if
!
!AA406 start
!AA601 sub
 if ((Domain%tipo == "bsph").and.(DBSPH%n_w > 0)) then
     write (nout,fmt101) &
  "  Wall velocity|",minvelo_w,"|",minlocvelo_w,"|",pg_w(minlocvelo_w)%coord(1),"|",pg_w(minlocvelo_w)%coord(2),"|",pg_w(minlocvelo_w)%coord(3), &
  "||",maxvelo_w,"|",maxlocvelo_w,"|",pg_w(maxlocvelo_w)%coord(1),"|",pg_w(maxlocvelo_w)%coord(2),"|",pg_w(maxlocvelo_w)%coord(3),"|"
     write (nout,fmt102) &
  "  Wall pressure|",minpres_w,"|",minlocpres_w,"|",pg_w(minlocpres_w)%coord(1),"|",pg_w(minlocpres_w)%coord(2),"|",pg_w(minlocpres_w)%coord(3), &
  "||",maxpres_w,"|",maxlocpres_w,"|",pg_w(maxlocpres_w)%coord(1),"|",pg_w(maxlocpres_w)%coord(2),"|",pg_w(maxlocpres_w)%coord(3),"|"
  endif
!AA406 end
!

!AA501b start
  if (n_bodies > 0) then
     write (nout,fmt101) &
  "Body part.vel. |",minvelo_bp,"|",minlocvelo_bp,"|",bp_arr(minlocvelo_bp)%pos(1),"|",bp_arr(minlocvelo_bp)%pos(2),"|", &
  bp_arr(minlocvelo_bp)%pos(3),"||",maxvelo_bp,"|",maxlocvelo_bp,"|",bp_arr(maxlocvelo_bp)%pos(1),"|", &
  bp_arr(maxlocvelo_bp)%pos(2),"|",bp_arr(maxlocvelo_bp)%pos(3),"|"
     write (nout,fmt102) &
  "Body part.pres.|",minpres_bp,"|",minlocpres_bp,"|",bp_arr(minlocpres_bp)%pos(1),"|",bp_arr(minlocpres_bp)%pos(2),"|", &
  bp_arr(minlocpres_bp)%pos(3),"||",maxpres_bp,"|",maxlocpres_bp,"|",bp_arr(maxlocpres_bp)%pos(1),"|", &
  bp_arr(maxlocpres_bp)%pos(2),"|",bp_arr(maxlocpres_bp)%pos(3),"|"
      write (nout,fmt101) &
  "Body velocity  |",minvelo_body,"|",minlocvelo_body,"|",body_arr(minlocvelo_body)%x_CM(1),"|", &
  body_arr(minlocvelo_body)%x_CM(2),"|",body_arr(minlocvelo_body)%x_CM(3),"||",maxvelo_body,"|",maxlocvelo_body,"|", &
  body_arr(maxlocvelo_body)%x_CM(1),"|",body_arr(maxlocvelo_body)%x_CM(2),"|",body_arr(maxlocvelo_body)%x_CM(3),"|"
      write (nout,fmt101) &
  "Body omega     |",minomega_body,"|",minlocomega_body,"|",body_arr(minlocomega_body)%x_CM(1),"|", &
  body_arr(minlocomega_body)%x_CM(2),"|",body_arr(minlocomega_body)%x_CM(3),"||",maxomega_body,"|",maxlocomega_body,"|", &
  body_arr(maxlocomega_body)%x_CM(1),"|",body_arr(maxlocomega_body)%x_CM(2),"|",body_arr(maxlocomega_body)%x_CM(3),"|"
  endif
!AA501b end

!AA504 start
  if ((Granular_flows_options%erosion_flag.ne.1).and.(Granular_flows_options%ID_erosion_criterion==1)) then
     write (nout,fmt101) &
  "tau_tauc       |",mintau_tauc,"|",minloctau_tauc,"|",pg(minloctau_tauc)%coord(1),"|",pg(minloctau_tauc)%coord(2),"|", &
  pg(minloctau_tauc)%coord(3),"||",maxtau_tauc,"|",maxloctau_tauc,"|",pg(maxloctau_tauc)%coord(1),"|", &
  pg(maxloctau_tauc)%coord(2),"|",pg(maxloctau_tauc)%coord(3),"|"
     write (nout,fmt101) &
  "u_star         |",minu_star,"|",minlocu_star,"|",pg(minlocu_star)%coord(1),"|",pg(minlocu_star)%coord(2),"|", &
  pg(minlocu_star)%coord(3),"||",maxu_star,"|",maxlocu_star,"|",pg(maxlocu_star)%coord(1),"|", &
  pg(maxlocu_star)%coord(2),"|",pg(maxlocu_star)%coord(3),"|"  
  if (Granular_flows_options%ID_erosion_criterion == 1) then       
     write (nout,fmt101) &
  "k_BetaGamma    |",mink_BetaGamma,"|",minlock_BetaGamma,"|",pg(minlock_BetaGamma)%coord(1),"|",pg(minlock_BetaGamma)%coord(2),"|", &
  pg(minlock_BetaGamma)%coord(3),"||",maxk_BetaGamma,"|",maxlock_BetaGamma,"|",pg(maxlock_BetaGamma)%coord(1),"|", &
  pg(maxlock_BetaGamma)%coord(2),"|",pg(maxlock_BetaGamma)%coord(3),"|"
  endif   
  endif
!AA504 end

!AA405 start
  else
     write (nout,'(128("."))') 
     write (nout,'(a)') "No particles inside the domain at the moment"
     write (nout,'(128("."))') 
  endif
!AA405 end
!
  it_print = it
!
  return
  end subroutine Print_Results
!---split

