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
! Program unit: Time_module            
! Description: Module for time recording.                                           
!-------------------------------------------------------------------------------
module Time_module
use Static_allocation_module
integer(4),parameter :: redcard = 20 ! Reduced length of strings
integer(4),parameter :: numb_subr = 21 ! Number of subroutine to count + 1
integer(4) :: ifiout
integer(4),dimension(numb_subr) :: tot_call
double precision, dimension(numb_subr,2) :: tot_times
character(LEN=40), dimension(numb_subr) :: tot_routines
character(len=redcard) :: routine
character(len=lencard) :: prefix,date_exec,nomefile_sta,                       &
                          nomecaso_rest,nomefile_inp,nomefile_msh,nomefile_cad,&
                          nomefile_sec,nomefile_spy,nomefile_res,nomefile_txt, &
                          nomefile_chk,nomefile_neu,case_title,case_subtitle
character(len=32) :: host_name
character(LEN=16) :: case_data  
character(LEN=8) :: case_hour
type iiar                                 
   sequence
   integer(4) ihr
   integer(4) imin
   integer(4) isec
end type
type(iiar) itime_struct
common /strings   / tot_routines,routine,date_exec,prefix,nomefile_sta,        &
                    nomecaso_rest,nomefile_inp,nomefile_msh,nomefile_cad,      &
                    nomefile_sec,nomefile_spy,host_name,nomefile_res,          &
                    nomefile_txt,nomefile_chk,nomefile_neu,case_data,case_hour,&
                    case_title,case_subtitle
common /timetable / tot_times,tot_call
data tot_routines(1) / "Initialization                          " /
data tot_routines(2) / "Gest_Input                              " /
data tot_routines(3) / "Gest_Trans                              " /
data tot_routines(4) / "   Loop_ghost                           " /
data tot_routines(5) / "   Loop_Irre_2D/3D                      " /
data tot_routines(6) / "     Motion Equation & Update velocity  " /
data tot_routines(7) / "     Velocity smoothing                 " /
data tot_routines(8) / "     Update position                    " /
data tot_routines(9) / "     Outgone/Source/Ordgrid             " /
data tot_routines(10) / "     Search for neighbourhood particles " /
data tot_routines(11) / "     Boundary integrals                 " /
data tot_routines(12) / "     Continuity equation                " /
data tot_routines(13) / "     State equation                     " /
data tot_routines(14) / "     Pressure smoothing                 " /
data tot_routines(15) / "     Apparent viscosity                 " /
data tot_routines(16) / "     Diffusion model                    " /
data tot_routines(17) / "     RK time integration                " /
data tot_routines(18) / "     Wall parameter update              " /
data tot_routines(19) / "     Rigid body transport               " /
data tot_routines(20) / "     Mixture viscosity                  " /
data tot_routines(numb_subr) / "Totale                                  " /
end module Time_module
 
