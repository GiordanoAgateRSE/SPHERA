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
! Program unit: SA_SPH_module            
! Description: Module for the semi-analytic approach (boundary treatment scheme) of Di Monaco et al. (2011, EACFM).                     
!----------------------------------------------------------------------------------------------------------------------------------

module SA_SPH_module
use Static_allocation_module
integer(4) :: iBIT,jBIT
double precision,public,parameter :: DELTAX_BIT = 0.1d0 
double precision,dimension(1:22,0:5) :: BoundIntegralTab2D
data ((BoundIntegralTab2D(iBIT,jBIT),jBIT=0,5),iBIT=1,22)  /  &
! delta_i     Int_Wd0s0    Int_IWd1s0   Int_Wd1s0    Int_Wd3s0    Int_Wd1s2
-9.900000E-02,6.916001E-01,5.671769E-01,1.614238E+00,7.329450E-01,8.813822E-01,&
 1.000000E-03,6.820633E-01,4.993001E-01,1.426603E+00,7.142813E-01,7.123606E-01,&
 1.010000E-01,6.725266E-01,4.314232E-01,1.238968E+00,6.956176E-01,5.433391E-01,&
 2.010000E-01,6.452481E-01,3.653947E-01,1.057834E+00,6.492812E-01,4.085521E-01,&
 3.010000E-01,6.025034E-01,3.028878E-01,8.871524E-01,5.851007E-01,3.020341E-01,&
 4.010000E-01,5.473293E-01,2.453079E-01,7.296175E-01,5.106232E-01,2.189841E-01,&
 5.010000E-01,4.832434E-01,1.937172E-01,5.872117E-01,4.319031E-01,1.553082E-01,&
 6.010000E-01,4.139976E-01,1.488267E-01,4.612826E-01,3.538107E-01,1.074685E-01,&
 7.010000E-01,3.434138E-01,1.109601E-01,3.526284E-01,2.802530E-01,7.237218E-02,&
 8.010000E-01,2.751568E-01,8.006434E-02,2.615428E-01,2.142205E-01,4.732315E-02,&
 9.010000E-01,2.125249E-01,5.573915E-02,1.877770E-01,1.577964E-01,2.998003E-02,&
 1.001000E+00,1.581202E-01,3.728314E-02,1.304301E-01,1.120873E-01,1.834320E-02,&
 1.101000E+00,1.131765E-01,2.379696E-02,8.756287E-02,7.679627E-02,1.076657E-02,&
 1.201000E+00,7.737532E-02,1.434262E-02,5.638309E-02,5.038764E-02,5.995379E-03,&
 1.301000E+00,4.995510E-02,8.042085E-03,3.441001E-02,3.129118E-02,3.118804E-03,&
 1.401000E+00,2.994721E-02,4.104289E-03,1.956203E-02,1.807987E-02,1.482120E-03,&
 1.501000E+00,1.623613E-02,1.842642E-03,1.008625E-02,9.464607E-03,6.216584E-04,&
 1.601000E+00,7.615461E-03,6.873181E-04,4.510624E-03,4.293336E-03,2.172841E-04,&
 1.701000E+00,2.842820E-03,1.913104E-04,1.609094E-03,1.552266E-03,5.682598E-05,&
 1.801000E+00,6.998383E-04,3.120105E-05,3.793334E-04,3.706037E-04,8.729655E-06,&
 1.901000E+00,6.214235E-05,1.372448E-06,3.231651E-05,3.195398E-05,3.625237E-07,&
 2.001000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00,0.000000E+00/
end module
    
