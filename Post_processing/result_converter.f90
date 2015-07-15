!cfile result_converter.f90
!************************************************************************************
!                             S P H E R A 6.0.0 
!
!                      Smoothed Particle Hydrodynamics Code
!
!************************************************************************************
!
! File name     : result_converter
!
! Last updating : November 23, 2011
!
! Improvement traceback:
!
! 00  R. Guandalini  31/07/2007   Module creation
! 01  G. Agate       24/10/2007   Inserted in SPHERA
! 02  Amicarelli     23/11/2011   multiple inlet
! 03  Amicarelli/Agate 30Nov11    BSPH: wall element parameters
!AA501b comment
! 04  Amicarelli-Agate 13nov12    Body dynamics
!AA504
! 05  Amicarelli       08Apr14    (v5.04) Modifications for granular flows and 3D erosion criterion
!AA601
! 06  Amicarelli       26Jan15    DBSPH-input (AA601). New DBSPH PV output. 
!
!************************************************************************************
! Module purpose : Converts the results into VTK distributions 
!
! Calling modules: Loop_Irre_2D, Loop_Irre_3D, diagnostic
!
! Called modules : none
!      
!************************************************************************************

  subroutine result_converter (str)
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
!.. Local Parameters
!  double precision, parameter    :: tol = 1.0d-04             ! tollerance
!
!.. Formal Arguments ..
  character(6),intent(IN) :: str
!
!.. Local Scalars ..
  character(len=256) :: stringa,header_string
  character(len=120) :: filevtk,prefix
  character(len=10)  :: cargo
  integer(4)         :: npi,i,j,k,k1,k2
  integer(4)         :: numcells,numpoints
  double precision   :: curtime
!
!.. Local Arrays ..
!  double precision, dimension(16)       :: appo
  integer(4), dimension(24)             :: iappo
  integer(4), dimension(:), allocatable :: finger

!.. executable statements
!
!.. check for time sampling is active
!
!!!    curtime = tempo
    curtime = tempo - MOD(tempo,abs(freq_time))  ! arrotondamento tempo corrente per animazioni
    if (curtime < val_time .and. index(str,'fine')==0) return
!    if (str == 'loop__') then
!      if (abs(val_time - curtime) > tol .and. abs(freq_time) > zero)  return
!      if (abs(val_time - curtime) > tol .and. freq_time < zero ) return
!    end if
!AA405
    if (nag > 0) then
!
    block = block + 1
    nblocchi = nblocchi + 1
    if (nblocchi > maxnumblock) then
      write (nscr,'(a)') ' ATTENZIONE !! nblocchi > maxnumblock in routine result_converter per file VTK.'
      write (nscr,'(a)') '               aumentare maxnumblock oppure diminuire frequenza memorizzazione per file VTK.'
      write (nout,'(a)') ' ATTENZIONE !! nblocchi > maxnumblock in routine result_converter per file VTK.'
      write (nout,'(a)') '               aumentare maxnumblock oppure diminuire frequenza memorizzazione per file VTK.'
      nblocchi = maxnumblock
    end if
    blocchi(nblocchi) = block
    Time_Block(nblocchi) = curtime
    prefix = nomecaso
!
!.. open the VTK formatted file VTKConverter_<casename>_<block>.vtk for the results storing
!
!!    write (cargo,'(f10.4)') curtime
!!    cargo = adjustl(cargo)
!!    i = scan(cargo,'.')
!!    cargo(i:(len_trim(cargo)-1)) = cargo(i+1:(len_trim(cargo)))
!!    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_time_"//cargo(1:len_trim(cargo))//".vtk"
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
!
!.. initialize the VTK formatted file
!
!.. cells are: 1 = vtk_poly_vertex including all the particle triplets in xyz cartesian system
!
    numcells = 1
!
!.. write the heading records
!
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
!
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
!
!.. evaluates and write the point coordinates 
!
    numpoints = count(pg(1:nag)%cella > 0)
    allocate (finger(numpoints))
    k = 0
    do npi = 1,nag
      if (pg(npi)%cella == 0) cycle
      k = k + 1
      finger(k) = npi
    end do
!
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
!
!.. write the coordinates of all the active particles
!
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (pg(finger(k))%coord(1),pg(finger(k))%coord(2),pg(finger(k))%coord(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
!
!.. writes the topology 
!
!.. write the first cell: vtk_poly_vertex of type 2
!
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the cell offsets
!
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the cell types
!
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
!
    write(unitvtk,'(a)') '    </Cells>'
!
!.. write the results on the VTK format file for post processing as point data
!
    write(unitvtk,'(a)') '<PointData>'
!
!AA501 rm start
!.. write the total velocity distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Total velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          j = 0
!          do k = k1,k2
!            j = j + 1
!            appo(j) = Dsqrt(pg(finger(k))%vel(1)*pg(finger(k))%vel(1) + pg(finger(k))%vel(2)*pg(finger(k))%vel(2) + &
!                      pg(finger(k))%vel(3)*pg(finger(k))%vel(3))
!          end do
!          write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!.. write the total velocity as vector distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
        do i = 1, numpoints,6
          k1 = i
          k2 = k1 + 5
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%vel(1),pg(finger(k))%vel(2),pg(finger(k))%vel(3),k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!.. write the X velocity component distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="X velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%vel(1),k = k1,k2)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the Y velocity component distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Y velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%vel(2),k = k1,k2)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the Z velocity component distribution
!
!      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Z velocity (m/s)" format="ascii" >'
!        do i = 1, numpoints,16
!          k1 = i
!          k2 = k1 + 15
!          if (k2 > numpoints) k2 = numpoints
!          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%vel(3),k = k1,k2)
!        end do
!      write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!.. write the pressure distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%pres,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the density distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Density (kg/mc)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%dens,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the mass distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Mass (kg)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%mass,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the viscosity distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Viscosity (mq/s)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%visc,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the Volume Fraction distribution
!
      if (diffusione) then
        write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Volume Fraction (mc/mc)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%VolFra,k = k1,k2)
        end do
        write(unitvtk,'(a)') '      </DataArray>'
      end if
!
!.. write the Internal Energy distribution
!
      if (esplosione) then
        write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Internal Energy (J/kg)" format="ascii" >'
        do i = 1, numpoints,16
          k1 = i
          k2 = k1 + 15
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%IntEn,k = k1,k2)
        end do
        write(unitvtk,'(a)') '      </DataArray>'
      end if
!
!.. write the medium distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Medium" format="ascii" >'
        do i = 1, numpoints,24
          k1 = i
          k2 = k1 + 23
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,24(1x,i6))') (pg(finger(k))%imed,k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the particle state distribution
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Status" format="ascii" >'
        do i = 1, numpoints,24
          k1 = i
          k2 = k1 + 23
          if (k2 > numpoints) k2 = numpoints
          iappo = 0
          j = 0
          do k = k1,k2
            j = j + 1
            if ( (pg(finger(k))%state == 'flu') .and. (index(Med(pg(finger(k))%imed)%tipo,"liquid") > 0) ) then
              iappo(k-k1+1) = 1
            else if ( (pg(finger(k))%state == 'flu') .and. (index(Med(pg(finger(k))%imed)%tipo,"granular") > 0) ) then
              iappo(k-k1+1) = 2
            else if ( (pg(finger(k))%state == 'sol') .and. (index(Med(pg(finger(k))%imed)%tipo,"granular") > 0) ) then
              iappo(k-k1+1) = 3
            else if ( (pg(finger(k))%state == 'flu') .and. (index(Med(pg(finger(k))%imed)%tipo,"gas") > 0) ) then
              iappo(k-k1+1) = 4
            end if
          end do
          write(unitvtk,'(8x,24(1x,i6))') (iappo(k),k = 1,j)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!.. write the particle position in the pg array
!
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Finger" format="ascii" >'
        do i = 1, numpoints,24
          k1 = i
          k2 = k1 + 23
          if (k2 > numpoints) k2 = numpoints
          write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
        end do
      write(unitvtk,'(a)') '      </DataArray>'
!
!AA406 start
 if (Domain%tipo == "bsph") then
!Writing Shepard's coefficient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Shepard coefficient" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%uni,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing the discrete Shepard's coefficient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="discrete Shepard" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%sigma,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing the integral Shepard's coefficient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="integral Shepard" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Gamma,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing the free surface flag
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="free surface" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,i8))') (pg(finger(k))%FS,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start

!
!Writing the absolute value of the density gradient
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Density gradient (kg/(m^4))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%drho(1)*pg(finger(k))%drho(1) + pg(finger(k))%drho(2)*pg(finger(k))%drho(2) + &
!                          pg(finger(k))%drho(3)*pg(finger(k))%drho(3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the density gradient
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="drho vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%drho(1),pg(finger(k))%drho(2),pg(finger(k))%drho(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!
!Writing the absolute value of the gradient of the x-component of velocity (d(u_x))
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_x) (1/(s))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%dvel(1,1)*pg(finger(k))%dvel(1,1) + pg(finger(k))%dvel(1,2)*pg(finger(k))%dvel(1,2) + &
!                          pg(finger(k))%dvel(1,3)*pg(finger(k))%dvel(1,3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the gradient of the x-component of velocity
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_x) vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(1,1),pg(finger(k))%dvel(1,2),pg(finger(k))%dvel(1,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    if (ncord==3) then
!
!AA501 rm start
!
!Writing the absolute value of the gradient of the y-component of velocity (d(u_y))
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_y) (1/(s))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%dvel(2,1)*pg(finger(k))%dvel(2,1) + pg(finger(k))%dvel(2,2)*pg(finger(k))%dvel(2,2) + &
!                          pg(finger(k))%dvel(2,3)*pg(finger(k))%dvel(2,3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the of the gradient of the y-component of velocity
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_y) vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(2,1),pg(finger(k))%dvel(2,2),pg(finger(k))%dvel(2,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    endif
!
!AA501 rm start
!
!Writing the absolute value of the gradient of the z-component of velocity (d(u_z))
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_z) (1/(s))" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg(finger(k))%dvel(3,1)*pg(finger(k))%dvel(3,1) + pg(finger(k))%dvel(3,2)*pg(finger(k))%dvel(3,2) + &
!                          pg(finger(k))%dvel(3,3)*pg(finger(k))%dvel(3,3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
!Writing the vector of the of the gradient of the z-component of velocity
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="d(u_z) vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(3,1),pg(finger(k))%dvel(3,2),pg(finger(k))%dvel(3,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    endif ! Domain%tipo == "bsph"
!AA406 end

!AA504 start
    if (Granular_flows_options%ID_erosion_criterion>=1) then
!Writing sigma_prime
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="sigma_prime" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%sigma_prime,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing sec_inv
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="sqrt_I2_eij" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%secinv,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing blt_flag
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="blt_flag" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,i8))') (pg(finger(k))%blt_flag,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'  
!Writing z-coordinate for top view 2D field of free surface
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="z(m)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%coord(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
!Writing Bingham number
    if (Granular_flows_options%viscosity_blt_formula==4) then
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Bn" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Bn,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    endif
    if (Granular_flows_options%erosion_flag.ne.1) then    
!Writing Beta
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Beta(radians)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Beta_slope,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
!Writing the vector normal_int
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="normal_int vector"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%normal_int(1),pg(finger(k))%normal_int(2),pg(finger(k))%normal_int(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing tau_tauc
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="tau_tauc" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%tau_tauc,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing u_star
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="u_star(m/s)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%u_star,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'  
    if (Granular_flows_options%ID_erosion_criterion == 1) then  
    if (ncord==3) then    
!Writing C_D
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="C_D" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%C_D,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
!Writing C_L
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="C_L" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%C_L,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!Writing Gamma
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Gamma(radians)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Gamma_slope,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>' 
    endif
!Writing k_BetaGamma
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="k_BetaGamma" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%k_BetaGamma,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'  
    endif
    endif
    endif 
!AA504 end

    write(unitvtk,'(a)') '    </PointData>'
!.. write the distribution data on cells
!
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
!
!.. close the result VTK file
!
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
!
!.. Flush of file contents is forced
!
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)
!
!AA406 rm !!! (moved later on)
!AA405
!    endif  !nag>0
!
!AA406 start
    if ((DBSPH%n_w > 0) .and. (Domain%tipo == "bsph")) then
!
! Open the VTK unstructured grid formatted file VTKConverter_<casename>_wall_<block>.vtk for the results storing
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_wall_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
! Initialization of the VTK formatted file (1 cell formally groupes all the points)
! Writing the heading records
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Evaluating and writing the point coordinates 
    numpoints = count(pg_w(1:DBSPH%n_w)%cella > 0)
    allocate (finger(numpoints))
    k = 0
    do npi = 1,DBSPH%n_w
      if (pg_w(npi)%cella == 0) cycle
      k = k + 1
      finger(k) = npi

    end do
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Writing the coordinates of all the particles
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (pg_w(finger(k))%coord(1),pg_w(finger(k))%coord(2),pg_w(finger(k))%coord(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
! Writing the topology 
! Writing the first cell: vtk_poly_vertex of type 2
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell offset
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell type
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Cells>'
! Writing the results on the VTK format file for post processing as point data
    write(unitvtk,'(a)') '<PointData>'
!
!AA501 rm start
!
! Writing the absolute value of velocity
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity (m/s)" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg_w(finger(k))%vel(1)*pg_w(finger(k))%vel(1) + pg_w(finger(k))%vel(2)*pg_w(finger(k))%vel(2) + &
!                    pg_w(finger(k))%vel(3)*pg_w(finger(k))%vel(3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
! Writing the velocity vector 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg_w(finger(k))%vel(1),pg_w(finger(k))%vel(2),pg_w(finger(k))%vel(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!
! Writing the absolute value of the normal (it should be the unity)
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="normal norm" format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       j = 0
!       do k = k1,k2
!          j = j + 1
!          appo(j) = Dsqrt(pg_w(finger(k))%normal(1)*pg_w(finger(k))%normal(1) + pg_w(finger(k))%normal(2)*pg_w(finger(k))%normal(2) + &
!                    pg_w(finger(k))%normal(3)*pg_w(finger(k))%normal(3))
!       end do
!       write(unitvtk,'(8x,16(1x,e12.5))') (appo(k),k = 1,j)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
! Writing the normal vectors n
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Normal vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg_w(finger(k))%normal(1),pg_w(finger(k))%normal(2),pg_w(finger(k))%normal(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm start
!
! Writing n_x
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="n_x " format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%normal(1),k = k1,k2)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
! Writing n_y
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="n_y " format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!       k2 = k1 + 15
!       if (k2 > numpoints) k2 = numpoints
!       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%normal(2),k = k1,k2)
!    end do
!    write(unitvtk,'(a)') '      </DataArray>'
! Writing n_z
!    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="n_z " format="ascii" >'
!    do i = 1, numpoints,16
!       k1 = i
!      k2 = k1 + 15
!      if (k2 > numpoints) k2 = numpoints
!      write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%normal(3),k = k1,k2)
!   end do
!    write(unitvtk,'(a)') '      </DataArray>'
!
!AA501 rm end
!
! Writing pressure
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%pres,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing weight
    if (ncord == 2) then
       write(unitvtk,'(a)') '      <DataArray type="Float32" Name="weight (m)" format="ascii" >'
       else
          write(unitvtk,'(a)') '      <DataArray type="Float32" Name="weight (m^2)" format="ascii" >'
    endif
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%weight,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing semi-particle mass 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="semi-particle mass (kg)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%mass,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
!AA601 start
! Writing semi-particle k_d 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="semi-particle k_d" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%k_d,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="semi-particle volume (m^3)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%volume,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'    
!AA601 end
! Write the particle position in the pg_w array
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Finger" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </PointData>'
! Writing the distribution data on cells
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
! Closing the VTK file
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly (FORTRAN 2003)
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)
!   
    endif  ! DBSPH%n_w>0 and Domain%tipo == "bsph"
!

!AA501b start
!Body dynamics option
    if (n_bodies > 0) then

! Body particles
! Open the VTK unstructured grid formatted file VTKConverter_<casename>_body-part_<block>.vtk for the results storing
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body-part_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
! Initialization of the VTK formatted file (1 cell formally groupes all the points)
! Writing the heading records
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Evaluating and writing the point coordinates 
    numpoints = count(bp_arr(1:n_body_part)%cell > 0)
    allocate (finger(numpoints))
    k = 0
    do npi = 1,n_body_part
      if (bp_arr(npi)%cell == 0) cycle
      k = k + 1
      finger(k) = npi
    end do
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Writing the coordinates of all the particles
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (bp_arr(finger(k))%pos(1),bp_arr(finger(k))%pos(2),bp_arr(finger(k))%pos(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
! Writing the topology 
! Writing the first cell: vtk_poly_vertex of type 2
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell offset
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell type
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Cells>'
! Writing the results on the VTK format file for post processing as point data
    write(unitvtk,'(a)') '<PointData>'
! Writing the velocity vector 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (bp_arr(finger(k))%vel(1),bp_arr(finger(k))%vel(2),bp_arr(finger(k))%vel(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the normal vectors n
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Normal vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (bp_arr(finger(k))%normal(1),bp_arr(finger(k))%normal(2), &
                                             bp_arr(finger(k))%normal(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing pressure
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%pres,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing mass 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="mass (kg)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%mass,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing area 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="area(m^2)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%area,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Write the particle position in the bp_arr array
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="Finger" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Write the body identifier 
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="Body" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (bp_arr(finger(k))%body,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </PointData>'
! Writing the distribution data on cells
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
! Closing the VTK file
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly (FORTRAN 2003)
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)

! Bodies
! Open the VTK unstructured grid formatted file VTKConverter_<casename>_body_<block>.vtk for the results storing
    write (cargo,'(i6)') block
    cargo = adjustl(cargo)
    filevtk = "VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body_"//cargo(1:len_trim(cargo))//".vtu"
    write (nscr,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    write (nout,'(a)') "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
    open (unit=unitvtk,file=filevtk,form='formatted',access='sequential',status='unknown')
    rewind (unit=unitvtk)
! Initialization of the VTK formatted file (1 cell formally groupes all the points)
! Writing the heading records
    write(unitvtk,'(a)') '<?xml version="1.0"?>'
    write(unitvtk,'(a)') &
    '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
    write(unitvtk,'(a)') '<UnstructuredGrid>'
    write(cargo,'(i6)') int(curtime)
    cargo = adjustl(cargo)
    header_string = "case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Evaluating and writing the point coordinates 
    numpoints = n_bodies
    allocate (finger(numpoints))
    k = 0
    do npi = 1,n_bodies
      k = k + 1
      finger(k) = npi
    end do
    write(cargo,'(i6)') numpoints
    cargo = adjustl(trim(cargo))
    stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
    write(cargo,'(i6)') numcells
    cargo = adjustl(trim(cargo))
    stringa = stringa(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
    write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Writing the coordinates of all the particles
    write(unitvtk,'(a)') '    <Points>'
    write(unitvtk,'(a)') '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
    do i = 1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))') (body_arr(finger(k))%x_CM(1),body_arr(finger(k))%x_CM(2), &
                                          body_arr(finger(k))%x_CM(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Points>'
! Writing the topology 
! Writing the first cell: vtk_poly_vertex of type 2
    write(unitvtk,'(a)') '    <Cells>'
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
    do i = 0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2 > numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell offset
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
    write(stringa,'(i8)') numpoints
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the cell type
    stringa = ' '
    write(unitvtk,'(a)') '      <DataArray type="UInt8" Name="types" format="ascii" >'
    stringa(1:6) = '     2'
    write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </Cells>'
! Writing the results on the VTK format file for post processing as point data
    write(unitvtk,'(a)') '<PointData>'
! Writing the velocity vector 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%u_CM(1),body_arr(finger(k))%u_CM(2), &
                                             body_arr(finger(k))%u_CM(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the angular velocity 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Angular velocity"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%omega(1),body_arr(finger(k))%omega(2), &
                                             body_arr(finger(k))%omega(3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing the moment of inertia
! first row 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Ic(x,:)"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(1,1),body_arr(finger(k))%Ic(1,2), &
                                             body_arr(finger(k))%Ic(1,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! second row 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Ic(y,:)"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(2,1),body_arr(finger(k))%Ic(2,2), &
                                             body_arr(finger(k))%Ic(2,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! third row 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="Ic(z,:)"  NumberOfComponents="3"  format="ascii" >'
    do i = 1, numpoints,6
       k1 = i
       k2 = k1 + 5
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(3,1),body_arr(finger(k))%Ic(3,2), &
                                             body_arr(finger(k))%Ic(3,3),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Writing mass 
    write(unitvtk,'(a)') '      <DataArray type="Float32" Name="mass (kg)" format="ascii" >'
    do i = 1, numpoints,16
       k1 = i
       k2 = k1 + 15
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,16(1x,e12.5))') (body_arr(finger(k))%mass,k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
! Write the particle position in the body_arr array
    write(unitvtk,'(a)') '      <DataArray type="Int32" Name="Finger" format="ascii" >'
    do i = 1, numpoints,24
       k1 = i
       k2 = k1 + 23
       if (k2 > numpoints) k2 = numpoints
       write(unitvtk,'(8x,24(1x,i6))') (finger(k),k = k1,k2)
    end do
    write(unitvtk,'(a)') '      </DataArray>'
    write(unitvtk,'(a)') '    </PointData>'
! Writing the distribution data on cells
    write(unitvtk,'(a)') '    <CellData>'
    write(unitvtk,'(a)') '    </CellData>'
! Closing the VTK file
    write(unitvtk,'(a)') '  </Piece>'
    write(unitvtk,'(a)') '</UnstructuredGrid>'
    write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly (FORTRAN 2003)
    flush(unitvtk)
    close (unitvtk)
    deallocate (finger)
 
    endif  
!AA501b end 
    
    endif  !nag>0
!AA406 end    

!.. increase the time sampling if active
!
     if (freq_time > zero) then
       val_time = val_time + abs(freq_time)
!
!AA404 
! truncation for compatibility with "curtime"
       val_time = val_time - MOD(val_time,abs(freq_time))
       if (val_time==curtime) val_time = val_time+abs(freq_time)
!
     else if (curtime >= val_time .and. val_time > zero) then
       val_time = 1.e20
     end if
!
!!  end do listing_loop
!
  return
  end subroutine result_converter
!---split

