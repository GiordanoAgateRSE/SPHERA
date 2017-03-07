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
! Program unit: result_converter               
! Description: Post-processing for .vtu (fluid dynamics parameters) and .vtk 
!              (geometry) files for Paraview. Subroutine call for the 
!              concatenation of the .txt temporary output files.      
!-------------------------------------------------------------------------------
subroutine result_converter(str)
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
character(6),intent(IN) :: str
integer(4) :: npi,i,j,k,k1,k2,numcells,numpoints
double precision :: curtime
integer(4),dimension(24) :: iappo
integer(4),dimension(:),allocatable :: finger
character(len=256) :: stringa,header_string
character(len=120) :: filevtk,prefix
character(len=10)  :: cargo
!------------------------
! Explicit interfaces
!------------------------
!------------------------
! Allocations
!------------------------
!------------------------
! Initializations
!------------------------
! .. check for time sampling is active
curtime = simulation_time  
!------------------------
! Statements
!------------------------
if ((curtime<val_time).and.(index(str,'fine')==0)) return
! To concatenate the ".txt" output files and remove the original ones
call cat_post_proc
! Rounding current physical time for Paraview .vtu files
   curtime = simulation_time - MOD(simulation_time,abs(freq_time))  
   if (nag>0) then
      block = block + 1
      nblocchi = nblocchi + 1
      if (nblocchi>maxnumblock) then
         write(nscr,'(a)')                                                     &
' Warning! nblocchi>maxnumblock in subroutine result_converter.'
         write(nscr,'(a)')                                                     &
         '    Increase maxnumblock or decrease output frequency for vtu files.'
      nblocchi = maxnumblock
   endif
   blocchi(nblocchi) = block
   Time_Block(nblocchi) = curtime
   prefix = nomecaso
   write(cargo,'(i8)') block
   cargo = adjustl(cargo)
   filevtk =                                                                   &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_"//cargo(1:len_trim(cargo))//".vtu"
   write(nout,'(a)')                                                           &
      "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
   open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',        &
      status='unknown')
   rewind(unit=unitvtk)
! Initializing the VTK formatted file
   numcells = 1
! Headings
   write(unitvtk,'(a)') '<?xml version="1.0"?>'
   write(unitvtk,'(a)')                                                        &
'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
   write(unitvtk,'(a)') '<UnstructuredGrid>'
   write(cargo,'(i8)') int(curtime)
   cargo = adjustl(cargo)
   header_string =                                                             &
"case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Coordinates 
   numpoints = count(pg(1:nag)%cella>0)
   allocate(finger(numpoints))
   k = 0
   do npi=1,nag
      if (pg(npi)%cella==0) cycle
      k = k + 1
      finger(k) = npi
   enddo
   write(cargo,'(i8)') numpoints
   cargo = adjustl(trim(cargo))
   stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
   write(cargo,'(i8)') numcells
   cargo = adjustl(trim(cargo))
   stringa = stringa                                                           &
(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
   write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Particle position 
   write(unitvtk,'(a)') '    <Points>'
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
   do i=1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2>numpoints) k2 = numpoints
      write (stringa,'(6(3(e12.5,1x)))')                                       &
         (pg(finger(k))%coord(1),pg(finger(k))%coord(2),pg(finger(k))%coord(3),k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
   write(unitvtk,'(a)') '    </Points>'
! Topology 
   write(unitvtk,'(a)') '    <Cells>'
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
   do i=0,numpoints-1,30
      k1 = i
      k2 = k1 + 29
      if (k2>numpoints) k2 = numpoints
      write (stringa,'(30i8)') (k,k=k1,k2)
      stringa = adjustl(trim(stringa))
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! Cell offsets
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Int32" Name="offsets" format="ascii" >'
   write(stringa,'(i8)') numpoints
   write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
   write(unitvtk,'(a)') '      </DataArray>'
! Cell types
   stringa = ' '
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="UInt8" Name="types" format="ascii" >'
   stringa(1:6) = '     2'
   write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
   write(unitvtk,'(a)') '      </DataArray>'
   write(unitvtk,'(a)') '    </Cells>'
   write(unitvtk,'(a)') '<PointData>'
! Velocity 
   write(unitvtk,'(a)')                                                        &
'      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
   do i=1,numpoints,6
      k1 = i
      k2 = k1 + 5
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%vel(1),              &
         pg(finger(k))%vel(2),pg(finger(k))%vel(3),k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! Pressure 
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
   do i=1,numpoints,16
      k1 = i
      k2 = k1 + 15
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%pres,k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! Density
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" Name="Density (kg/mc)" format="ascii" >'
   do i=1,numpoints,16
      k1 = i
      k2 = k1 + 15
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%dens,k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! Mass
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" Name="Mass (kg)" format="ascii" >'
   do i=1,numpoints,16
      k1 = i
      k2 = k1 + 15
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%mass,k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! Viscosity
   write(unitvtk,'(a)')                                                        &
'      <DataArray type="Float32" Name="Viscosity (m^2/s)" format="ascii" >'
   do i=1,numpoints,16
      k1 = i
      k2 = k1 + 15
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%visc,k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! z-coordinate
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" Name="z(m)" format="ascii" >'
   do i=1,numpoints,16
      k1 = i
      k2 = k1 + 15
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%coord(3),k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>' 
! Volume Fraction
   if (diffusione) then
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Volume Fraction (mc/mc)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%VolFra,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
   endif
! Specific internal energy
   if (esplosione) then
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Internal Energy (J/kg)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%IntEn,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
   endif
! Fluid ID
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" Name="Medium" format="ascii" >'
   do i=1,numpoints,24
      k1 = i
      k2 = k1 + 23
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,24(1x,i6))') (pg(finger(k))%imed,k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! Particle "status"
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" Name="Status" format="ascii" >'
   do i=1,numpoints,24
      k1 = i
      k2 = k1 + 23
      if (k2>numpoints) k2 = numpoints
      iappo = 0
      j = 0
      do k=k1,k2
         j = j + 1
         if ((pg(finger(k))%state=="flu").and.                                 &
            (index(Med(pg(finger(k))%imed)%tipo,"liquid")>0)) then
            iappo(k-k1+1) = 1
            elseif ((pg(finger(k))%state=="flu").and.                          &
               (index(Med(pg(finger(k))%imed)%tipo,"granular")>0) ) then
               iappo(k-k1+1) = 2
               elseif ((pg(finger(k))%state=="sol").and.                       &
                  (index(Med(pg(finger(k))%imed)%tipo,"granular")>0) ) then
                  iappo(k-k1+1) = 3
                  elseif ((pg(finger(k))%state=="flu").and.                    &
                     (index(Med(pg(finger(k))%imed)%tipo,"gas")>0) ) then
                     iappo(k-k1+1) = 4
         endif
      enddo
      write(unitvtk,'(8x,24(1x,i6))') (iappo(k),k = 1,j)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
! Current particle ID
   write(unitvtk,'(a)')                                                        &
      '      <DataArray type="Float32" Name="Finger" format="ascii" >'
   do i=1,numpoints,24
      k1 = i
      k2 = k1 + 23
      if (k2>numpoints) k2 = numpoints
      write(unitvtk,'(8x,24(1x,i8))') (finger(k),k=k1,k2)
   enddo
   write(unitvtk,'(a)') '      </DataArray>'
   if (Domain%tipo=="bsph") then
! Shepard coefficient
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Shepard coefficient" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%uni,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Discrete Shepard coefficient
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Float32" Name="discrete Shepard" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%sigma,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Discrete Shepard coefficient in the same fluid
      if (NMedium>1) then
         write(unitvtk,'(a)')                                                  &
'      <DataArray type="Float32" Name="discrete Shepard same phase" format="ascii" >'
         do i=1,numpoints,16
            k1 = i
            k2 = k1 + 15
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%sigma_same_fluid,&
               k=k1,k2)
         enddo
         write(unitvtk,'(a)') '      </DataArray>'
      endif      
! Integral Shepard coefficient
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="integral Shepard" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%Gamma,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Free surface flag
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="UInt8" Name="free surface" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,i8))') (pg(finger(k))%FS,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Density gradient
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="drho vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%drho(1),          &
            pg(finger(k))%drho(2),pg(finger(k))%drho(3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Gradient of the x-component of velocity
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="d(u_x) vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(1,1),        &
            pg(finger(k))%dvel(1,2),pg(finger(k))%dvel(1,3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      if (ncord==3) then
! Gradient of the y-component of velocity
         write(unitvtk,'(a)')                                                  &
'      <DataArray type="Float32" Name="d(u_y) vectors"  NumberOfComponents="3"  format="ascii" >'
         do i=1,numpoints,6
            k1 = i
            k2 = k1 + 5
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(2,1),     &
               pg(finger(k))%dvel(2,2),pg(finger(k))%dvel(2,3),k=k1,k2)
         enddo
         write(unitvtk,'(a)') '      </DataArray>'
      endif
! Gradient of the z-component of velocity
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="d(u_z) vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%dvel(3,1),        &
            pg(finger(k))%dvel(3,2),pg(finger(k))%dvel(3,3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
   endif
! Laminar_flag
   if (.not.((Domain%tipo=="bsph").and.(DBSPH%slip_ID<2))) then 
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="laminar_flag" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,i8))') (pg(finger(k))%laminar_flag,k = k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
   endif
   if (Granular_flows_options%ID_erosion_criterion>=1) then
! sigma_prime_m
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="sigma_prime_m (Pa)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))')                                    &
            (pg(finger(k))%sigma_prime_m,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! pres_fluid
      write(unitvtk,'(a)') '      <DataArray type="Float32" Name="pres_fluid (Pa)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))')                                    &
            (pg(finger(k))%pres_fluid,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! sec_inv
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Float32" Name="sqrt_I2_eij (s^-1)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%secinv,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! blt_flag
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="blt_flag" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,i8))') (pg(finger(k))%blt_flag,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'  
! normal_int_mixture_top
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="normal_int_mixture_top vector"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
! Zeroing the mixture top normal, for those particles, which do not represent 
! this interface
         do k=k1,k2
            if (pg(finger(k))%blt_flag/=2)                                     &
               pg(finger(k))%normal_int_mixture_top(:) = 0.d0
         enddo
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (                                &
            pg(finger(k))%normal_int_mixture_top(1),                           &
            pg(finger(k))%normal_int_mixture_top(2),                           &
            pg(finger(k))%normal_int_mixture_top(3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! normal_int_sat_top
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="normal_int_sat_top vector"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
! Zeroing the mixture top normal, for those particles, which do not represent 
! this interface
         do k=k1,k2
            if (pg(finger(k))%blt_flag/=5)                                     &
               pg(finger(k))%normal_int_sat_top(:) = 0.d0
         enddo
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (                                &
            pg(finger(k))%normal_int_sat_top(1),                               &
            pg(finger(k))%normal_int_sat_top(2),                               &
            pg(finger(k))%normal_int_sat_top(3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      if (Granular_flows_options%erosion_flag.ne.1) then
! Beta
         write(unitvtk,'(a)')                                                  &
         '      <DataArray type="Float32" Name="Beta(radians)" format="ascii" >'
         do i=1,numpoints,16
            k1 = i
            k2 = k1 + 15
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,16(1x,e12.5))')                                 &
               (pg(finger(k))%Beta_slope,k=k1,k2)
         enddo
         write(unitvtk,'(a)') '      </DataArray>' 
! tau_tauc
         write(unitvtk,'(a)')                                                  &
'      <DataArray type="Float32" Name="tau_tauc" format="ascii" >'
         do i=1,numpoints,16
            k1 = i
            k2 = k1 + 15
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%tau_tauc,k=k1,k2)
         enddo
         write(unitvtk,'(a)') '      </DataArray>'
! u_star
         write(unitvtk,'(a)')                                                  &
         '      <DataArray type="Float32" Name="u_star(m/s)" format="ascii" >'
         do i=1,numpoints,16
            k1 = i
            k2 = k1 + 15
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%u_star,k=k1,k2)
         enddo
         write(unitvtk,'(a)') '      </DataArray>'
! normal_int
         write(unitvtk,'(a)')                                                  &
'      <DataArray type="Float32" Name="normal_int vector"  NumberOfComponents="3"  format="ascii" >'
         do i=1,numpoints,6
            k1 = i
            k2 = k1 + 5
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg(finger(k))%normal_int(1), &
               pg(finger(k))%normal_int(2),pg(finger(k))%normal_int(3),k=k1,k2)
         enddo
         write(unitvtk,'(a)') '      </DataArray>'    
         if (Granular_flows_options%ID_erosion_criterion==1) then  
            if (ncord==3) then    
! C_D
               write(unitvtk,'(a)')                                            &
            '      <DataArray type="Float32" Name="C_D" format="ascii" >'
               do i=1,numpoints,16
                  k1 = i
                  k2 = k1 + 15
                  if (k2>numpoints) k2 = numpoints
                  write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%C_D,k=k1,k2)
               enddo
               write(unitvtk,'(a)') '      </DataArray>' 
! C_L
               write(unitvtk,'(a)')                                            &
                  '      <DataArray type="Float32" Name="C_L" format="ascii" >'
               do i=1,numpoints,16
                  k1 = i
                  k2 = k1 + 15
                  if (k2>numpoints) k2 = numpoints
                  write(unitvtk,'(8x,16(1x,e12.5))') (pg(finger(k))%C_L,k=k1,k2)
               enddo
               write(unitvtk,'(a)') '      </DataArray>'
! Gamma
               write(unitvtk,'(a)')                                            &
'      <DataArray type="Float32" Name="Gamma(radians)" format="ascii" >'
               do i=1,numpoints,16
                  k1 = i
                  k2 = k1 + 15
                  if (k2>numpoints) k2 = numpoints
                  write(unitvtk,'(8x,16(1x,e12.5))')                           &
                     (pg(finger(k))%Gamma_slope,k=k1,k2)
               enddo
               write(unitvtk,'(a)') '      </DataArray>' 
            endif
! k_BetaGamma
            write(unitvtk,'(a)')                                               &
         '      <DataArray type="Float32" Name="k_BetaGamma" format="ascii" >'
            do i=1,numpoints,16
               k1 = i
               k2 = k1 + 15
               if (k2>numpoints) k2 = numpoints
               write(unitvtk,'(8x,16(1x,e12.5))')                              &
                 (pg(finger(k))%k_BetaGamma,k=k1,k2)
            enddo
            write(unitvtk,'(a)') '      </DataArray>'  
         endif
      endif
   endif 
   write(unitvtk,'(a)') '    </PointData>'
   write(unitvtk,'(a)') '    <CellData>'
   write(unitvtk,'(a)') '    </CellData>'
! Close the .vtu file
   write(unitvtk,'(a)') '  </Piece>'
   write(unitvtk,'(a)') '</UnstructuredGrid>'
   write(unitvtk,'(a)') '</VTKFile>'
! Flush unit content
   flush(unitvtk)
   close (unitvtk)
   deallocate(finger)
   if ((DBSPH%n_w>0).and.(Domain%tipo=="bsph")) then
! Open a vtu file for DB-SPH wall and semi-particle parameters  
! VTKConverter_<casename>_wall_<block>.vtk 
      write (cargo,'(i8)') block
      cargo = adjustl(cargo)
      filevtk =                                                                &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_wall_"//cargo(1:len_trim(cargo))//".vtu"
      write(nout,'(a)')                                                        &
         "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
      open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',     &
         status='unknown')
      rewind(unit=unitvtk)
! Initialization of the .vtu file (1 cell formally groupes all the 
! points)
! Heading
      write(unitvtk,'(a)') '<?xml version="1.0"?>'
      write(unitvtk,'(a)')                                                     &
'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      write(unitvtk,'(a)') '<UnstructuredGrid>'
      write(cargo,'(i8)') int(curtime)
      cargo = adjustl(cargo)
      header_string =                                                          &
"case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Point coordinates 
      numpoints = count(pg_w(1:DBSPH%n_w)%cella>0)
      allocate(finger(numpoints))
      k = 0
      do npi=1,DBSPH%n_w
         if (pg_w(npi)%cella==0) cycle
         k = k + 1
         finger(k) = npi
      enddo
      write(cargo,'(i8)') numpoints
      cargo = adjustl(trim(cargo))
      stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
      write(cargo,'(i8)') numcells
      cargo = adjustl(trim(cargo))
      stringa = stringa                                                        &
(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Position
      write(unitvtk,'(a)') '    <Points>'
      write(unitvtk,'(a)')                                                     &
      '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write (stringa,'(6(3(e12.5,1x)))') (pg_w(finger(k))%coord(1),         &
         pg_w(finger(k))%coord(2),pg_w(finger(k))%coord(3),k=k1,k2)
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </Points>'
! Topology 
      write(unitvtk,'(a)') '    <Cells>'
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
      do i=0,numpoints-1,30
         k1 = i
         k2 = k1 + 29
         if (k2>numpoints) k2 = numpoints
         write (stringa,'(30i8)') (k,k=k1,k2)
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Cell offset
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="offsets" format="ascii" >'
      write(stringa,'(i8)') numpoints
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      write(unitvtk,'(a)') '      </DataArray>'
! Cell type
      stringa = ' '
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="UInt8" Name="types" format="ascii" >'
      stringa(1:6) = '     2'
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </Cells>'
! Actual results: start
      write(unitvtk,'(a)') '<PointData>'
! Velocity 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg_w(finger(k))%vel(1),         &
            pg_w(finger(k))%vel(2),pg_w(finger(k))%vel(3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Normal vectors
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Normal vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (pg_w(finger(k))%normal(1),      &
            pg_w(finger(k))%normal(2),pg_w(finger(k))%normal(3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Velocity gradient in VSL (projected along the wall element normal) times the
! shear viscosity
      if (DBSPH%slip_ID>0) then
         write(unitvtk,'(a)')                                                  &
'      <DataArray type="Float32" Name="Velocity gradient in VSL"  NumberOfComponents="3"  format="ascii" >'
         do i=1,numpoints,6
            k1 = i
            k2 = k1 + 5
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,6(3(1x,e12.5)))') (                             &
               pg_w(finger(k))%grad_vel_VSL_times_mu(1),                       &
               pg_w(finger(k))%grad_vel_VSL_times_mu(2),                       &
               pg_w(finger(k))%grad_vel_VSL_times_mu(3),k=k1,k2)
         enddo
         write(unitvtk,'(a)') '      </DataArray>'
      endif
! Pressure
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%pres,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Weight
      if (ncord==2) then
         write(unitvtk,'(a)')                                                  &
            '      <DataArray type="Float32" Name="weight (m)" format="ascii" >'
         else
            write(unitvtk,'(a)')                                               &
'      <DataArray type="Float32" Name="weight (m^2)" format="ascii" >'
      endif
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%weight,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Semi-particle mass 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="semi-particle mass (kg)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%mass,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Semi-particle k_d 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="semi-particle k_d" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%k_d,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Semi-particle volume      
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="semi-particle volume (m^3)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (pg_w(finger(k))%volume,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Semi-particle kinematic viscosity
      if (DBSPH%slip_ID>0) then
         write(unitvtk,'(a)')                                                  &
'      <DataArray type="Float32" Name="semi-particle kinematic viscosity (m^2/s)" format="ascii" >'
         do i=1,numpoints,16
            k1 = i
            k2 = k1 + 15
            if (k2>numpoints) k2 = numpoints
            write(unitvtk,'(8x,16(1x,e12.5))')                                 &
               (pg_w(finger(k))%kin_visc_semi_part,k=k1,k2)
         enddo      
         write(unitvtk,'(a)') '      </DataArray>'    
      endif
! Wall element ID 
      write(unitvtk,'(a)')                                                     &
      '      <DataArray type="Float32" Name="Finger" format="ascii" >'
      do i=1,numpoints,24
         k1 = i
         k2 = k1 + 23
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,24(1x,i6))') (finger(k),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </PointData>'
! Distribution data on cells
      write(unitvtk,'(a)') '    <CellData>'
      write(unitvtk,'(a)') '    </CellData>'
! Closing the .vtu file
      write(unitvtk,'(a)') '  </Piece>'
      write(unitvtk,'(a)') '</UnstructuredGrid>'
      write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly
      flush(unitvtk)
      close (unitvtk)
      deallocate(finger)
   endif
   if (n_bodies>0) then
! Body Transport post-processing for .vtu files: start
! Body particles
! Open the .vtu unstructured grid formatted file 
! VTKConverter_<casename>_body-part_<block>.vtk for the results storing
      write (cargo,'(i8)') block
      cargo = adjustl(cargo)
      filevtk =                                                                &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body-part_"//cargo(1:len_trim(cargo))//".vtu"
      write (nout,'(a)')                                                       &
          "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
      open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',     &
         status='unknown')
      rewind(unit=unitvtk)
! Initialization of the .vtu formatted file (1 cell formally groupes all the
! points)
! Headings
      write(unitvtk,'(a)') '<?xml version="1.0"?>'
      write(unitvtk,'(a)')                                                     &
'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      write(unitvtk,'(a)') '<UnstructuredGrid>'
      write(cargo,'(i8)') int(curtime)
      cargo = adjustl(cargo)
      header_string =                                                          &
"case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Point coordinates 
      numpoints = count(bp_arr(1:n_body_part)%cell>0)
      allocate(finger(numpoints))
      k = 0
      do npi=1,n_body_part
         if (bp_arr(npi)%cell==0) cycle
         k = k + 1
         finger(k) = npi
      enddo
      write(cargo,'(i8)') numpoints
      cargo = adjustl(trim(cargo))
      stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
      write(cargo,'(i8)') numcells
      cargo = adjustl(trim(cargo))
      stringa = stringa                                                        &
(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Position
      write(unitvtk,'(a)') '    <Points>'
      write(unitvtk,'(a)')                                                     &
      '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write (stringa,'(6(3(e12.5,1x)))') (bp_arr(finger(k))%pos(1),         &
            bp_arr(finger(k))%pos(2),bp_arr(finger(k))%pos(3),k=k1,k2)
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </Points>'
! Topology 
      write(unitvtk,'(a)') '    <Cells>'
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
      do i=0,numpoints-1,30
         k1 = i
         k2 = k1 + 29
         if (k2>numpoints) k2 = numpoints
         write (stringa,'(30i8)') (k,k=k1,k2)
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Cell offset
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="offsets" format="ascii" >'
      write(stringa,'(i8)') numpoints
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      write(unitvtk,'(a)') '      </DataArray>'
! Cell type
      stringa = ' '
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="UInt8" Name="types" format="ascii" >'
      stringa(1:6) = '     2'
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </Cells>'
! Actual results: start 
      write(unitvtk,'(a)') '<PointData>'
! Velocity 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (bp_arr(finger(k))%vel(1),       &
            bp_arr(finger(k))%vel(2),bp_arr(finger(k))%vel(3),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Normal vectors
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Normal vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (bp_arr(finger(k))%normal(1),    &
                                               bp_arr(finger(k))%normal(2),    &
                                               bp_arr(finger(k))%normal(3),    &
                                               k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Pressure
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Float32" Name="Pressure (Pa)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%pres,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Mass 
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Float32" Name="mass (kg)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%mass,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Area 
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Float32" Name="area(m^2)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (bp_arr(finger(k))%area,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Particle ID 
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="Finger" format="ascii" >'
      do i=1,numpoints,24
         k1 = i
         k2 = k1 + 23
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,24(1x,i6))') (finger(k),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Body ID  
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="Body" format="ascii" >'
      do i=1,numpoints,24
         k1 = i
         k2 = k1 + 23
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,24(1x,i6))') (bp_arr(finger(k))%body,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </PointData>'
! Distribution data on cells
      write(unitvtk,'(a)') '    <CellData>'
      write(unitvtk,'(a)') '    </CellData>'
! Closing the .vtu file
      write(unitvtk,'(a)') '  </Piece>'
      write(unitvtk,'(a)') '</UnstructuredGrid>'
      write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly 
      flush(unitvtk)
      close (unitvtk)
      deallocate(finger)
! Bodies: start
! Open the .vtu unstructured grid formatted file 
! VTKConverter_<casename>_body_<block>.vtk for the results storing
      write (cargo,'(i8)') block
      cargo = adjustl(cargo)
      filevtk =                                                                &
"VTKConverter_"//prefix(1:len_trim(prefix))//"_block_body_"//cargo(1:len_trim(cargo))//".vtu"
      write(nout,'(a)')                                                        &
         "VTK formatted converted file  : "//filevtk(1:len_trim(filevtk))
      open(unit=unitvtk,file=filevtk,form='formatted',access='sequential',     &
         status='unknown')
      rewind(unit=unitvtk)
! Initialization of the .vtu formatted file 
! (1 cell formally groupes all the points)
! Headings
      write(unitvtk,'(a)') '<?xml version="1.0"?>'
      write(unitvtk,'(a)')                                                     &
'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian" compressor="vtkZLibDataCompressor">'
      write(unitvtk,'(a)') '<UnstructuredGrid>'
      write(cargo,'(i8)') int(curtime)
      cargo = adjustl(cargo)
      header_string =                                                          &
"case "//prefix(1:len_trim(prefix))//" * time "//cargo(1:len_trim(cargo))//" (s)"
! Point coordinates 
      numpoints = n_bodies
      allocate(finger(numpoints))
      k = 0
      do npi=1,n_bodies
         k = k + 1
         finger(k) = npi
      enddo
      write(cargo,'(i8)') numpoints
      cargo = adjustl(trim(cargo))
      stringa = '  <Piece NumberOfPoints="'//cargo(1:len_trim(cargo))
      write(cargo,'(i8)') numcells
      cargo = adjustl(trim(cargo))
      stringa = stringa                                                        &
(1:len_trim(stringa))//'"        NumberOfCells="'//cargo(1:len_trim(cargo))//'"      >'
      write(unitvtk,'(a)') stringa(1:len_trim(stringa))
! Position
      write(unitvtk,'(a)') '    <Points>'
      write(unitvtk,'(a)')                                                     &
      '      <DataArray type="Float32" NumberOfComponents="3" format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write (stringa,'(6(3(e12.5,1x)))') (body_arr(finger(k))%x_CM(1),      &
                                             body_arr(finger(k))%x_CM(2),      &
                                             body_arr(finger(k))%x_CM(3),      &
                                             k=k1,k2)
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </Points>'
! Topology 
      write(unitvtk,'(a)') '    <Cells>'
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="connectivity" format="ascii" >'
      do i=0,numpoints-1,30
         k1 = i
         k2 = k1 + 29
         if (k2>numpoints) k2 = numpoints
         write (stringa,'(30i8)') (k,k=k1,k2)
         stringa = adjustl(trim(stringa))
         write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Cell offset
      write(unitvtk,'(a)') '      <DataArray type="Int32" Name="offsets" format="ascii" >'
      write(stringa,'(i8)') numpoints
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      write(unitvtk,'(a)') '      </DataArray>'
! Cell type
      stringa = ' '
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="UInt8" Name="types" format="ascii" >'
      stringa(1:6) = '     2'
      write(unitvtk,'(8x,a)') stringa(1:len_trim(stringa))
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </Cells>'
! Actual results: start 
      write(unitvtk,'(a)') '<PointData>'
! Velocity 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Velocity vectors"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%u_CM(1),    &
                                               body_arr(finger(k))%u_CM(2),    &
                                               body_arr(finger(k))%u_CM(3),    &
                                               k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Angular velocity 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Angular velocity"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%omega(1),   &
                                               body_arr(finger(k))%omega(2),   &
                                               body_arr(finger(k))%omega(3),   &
                                               k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Moment of inertia
! First row 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Ic(x,:)"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(1,1),    &
                                               body_arr(finger(k))%Ic(1,2),    &
                                               body_arr(finger(k))%Ic(1,3),    &
                                               k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Second row 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Ic(y,:)"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(2,1),    &
                                               body_arr(finger(k))%Ic(2,2),    &
                                               body_arr(finger(k))%Ic(2,3),    &
                                               k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Third row 
      write(unitvtk,'(a)')                                                     &
'      <DataArray type="Float32" Name="Ic(z,:)"  NumberOfComponents="3"  format="ascii" >'
      do i=1,numpoints,6
         k1 = i
         k2 = k1 + 5
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,6(3(1x,e12.5)))') (body_arr(finger(k))%Ic(3,1),    &
                                               body_arr(finger(k))%Ic(3,2),    &
                                               body_arr(finger(k))%Ic(3,3),    &
                                               k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Mass 
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Float32" Name="mass (kg)" format="ascii" >'
      do i=1,numpoints,16
         k1 = i
         k2 = k1 + 15
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,16(1x,e12.5))') (body_arr(finger(k))%mass,k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
! Body ID
      write(unitvtk,'(a)')                                                     &
         '      <DataArray type="Int32" Name="Finger" format="ascii" >'
      do i=1,numpoints,24
         k1 = i
         k2 = k1 + 23
         if (k2>numpoints) k2 = numpoints
         write(unitvtk,'(8x,24(1x,i6))') (finger(k),k=k1,k2)
      enddo
      write(unitvtk,'(a)') '      </DataArray>'
      write(unitvtk,'(a)') '    </PointData>'
! Distribution data on cells
      write(unitvtk,'(a)') '    <CellData>'
      write(unitvtk,'(a)') '    </CellData>'
! Closing the .vtu file
      write(unitvtk,'(a)') '  </Piece>'
      write(unitvtk,'(a)') '</UnstructuredGrid>'
      write(unitvtk,'(a)') '</VTKFile>'
! Flushing the unit explicitly 
      flush(unitvtk)
      close (unitvtk)
      deallocate(finger)
   endif  
endif 
! Updating the last output time for .vtu files 
if (freq_time>zero) then
   val_time = val_time + abs(freq_time)
   elseif ((curtime>=val_time).and.(val_time>zero)) then
      val_time = 1.e20
endif
!------------------------
! Deallocations
!------------------------
return
end subroutine result_converter

