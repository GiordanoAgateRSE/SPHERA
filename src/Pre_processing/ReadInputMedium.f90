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
! Program unit: ReadInputMedium                          
! Description: to read the input data from the section "medium" of SPHERA main 
!              input file.                       
!-------------------------------------------------------------------------------
subroutine ReadInputMedium(NumberEntities,Med,ainp,comment,nrighe,ier,ninp,    &
                           nout,nscr)
!------------------------
! Modules
!------------------------ 
use Static_allocation_module
use Hybrid_allocation_module
!------------------------
! Declarations
!------------------------
implicit none
integer(4) :: nrighe,ier,ninp,nout,nscr
integer(4),dimension(20) :: NumberEntities
type (TyMedium),dimension(NMedium) :: Med
character(1) :: comment
character(100) :: ainp
logical :: saturated_medium_flag
integer(4) :: index,nitersol,ioerr
double precision :: den0,eps,alfaMon,betaMon,visc,viscmx,taucri,cuin,phi,Cs
double precision :: cons,codif,Settling,coes,Rough,d50,Gamma,InitialIntEn,d_90
double precision :: limiting_viscosity
double precision :: porosity
character(8) :: tipo
character(100) :: token
character(100),external :: GetToken,lcase
logical,external :: ReadCheck
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
call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM DATA",ninp,nout)) return
do while (TRIM(lcase(ainp))/="##### end medium #####")
   read(ainp,*,iostat=ioerr) tipo
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM TYPE",ninp,nout)) return
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   read(ainp,*,iostat=ioerr) index
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM INDEX",ninp,nout)) return
   NumberEntities(2) = max(NumberEntities(2),index)
   den0 = zero
   eps = zero
   alfaMon = zero
   betaMon = zero
   codif = zero
   Settling = zero
   Gamma = zero
   InitialIntEn = zero
   visc = zero
   Rough = zero
   Cs = zero
   cons = zero
   viscmx = zero
   taucri = zero
   cuin = zero
   coes = zero
   limiting_viscosity = zero
   phi = zero
   saturated_medium_flag = .false.
   d50 = zero
   nitersol = 0
   porosity = 0.d0
   d_90 = 0.d0
   tipo = lcase(tipo)
   select case (tipo)
! Imposed dynamic viscosity for gases (explosion)
      case ("gas     ")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) den0,eps
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "WATER DENSITY & COMPRIMIBILITY",ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) alfaMon, betaMon
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",ninp,&
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) codif, Settling
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,nout))&
            return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) gamma, InitialIntEn
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,nout))&
            return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) visc
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DYNAMIC VISCOSITY",ninp,    &
            nout)) return            
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Rough
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,nout&
            )) return
      case ("liquid  ")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) den0,eps
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "WATER DENSITY & COMPRIMIBILITY",ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) alfaMon, betaMon
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",     &
            ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) codif,Settling
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) gamma,InitialIntEn
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) visc
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DYNAMIC VISCOSITY",ninp,    & 
            nout)) return            
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Rough
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,    &
            nout)) return
! Newtonian fluids with imposed dynamic turbulent viscosity (mixing-length 
! scheme, depending on dx, see "Smagorinsky") 
      case ("smagorin")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) den0,eps
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "WATER DENSITY & COMPRIMIBILITY",ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) alfaMon, betaMon
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",     &
            ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) codif, Settling
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) gamma,InitialIntEn
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) visc
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DYNAMIC VISCOSITY",ninp,    &
            nout)) return            
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Cs
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"SMAGORINSKY CONST",ninp,    &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Rough
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,    &
            nout)) return   
! Non-Newtoniani fluids with apparent viscosity (Chen formula)
      case ("general ")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) den0,eps
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "WATER DENSITY & COMPRIMIBILITY",ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) alfaMon, betaMon
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",     &
            ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) codif,Settling
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) gamma,InitialIntEn
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) cons,viscmx
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"CONSISTENCY & VISCO MAX",   &
            ninp,nout)) return            
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) taucri, cuin
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"TAUCRI & CUIN",ninp,nout)   &
            ) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) Rough
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ROUGH COEFFICIENT",ninp,    &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) nitersol
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"N° ITERATION SOLID",ninp,   &
            nout)) return
         visc = viscmx
! Non-Newtonian fluids with apparent viscosity 
      case ("granular")
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) den0,eps
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "WATER DENSITY & COMPRIMIBILITY",ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) alfaMon, betaMon
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"ALPHA e BETA MONAGHAN",     &
            ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) codif,Settling
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"DIFFUSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) gamma,InitialIntEn
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"EXPLOSION COEFF",ninp,      &
            nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) phi,saturated_medium_flag
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"PHI, SATURATED_MEDIUM_FLAG",&
            ninp,nout)) return
         call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
         read(ainp,*,iostat=ioerr) coes,viscmx,visc,limiting_viscosity
         if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                             &
            "COHESION, VISCO MAX, VISCO & LIMITING VISCOSITY",ninp,nout)) return            
         if (Granular_flows_options%ID_erosion_criterion==1) then
            call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
            read(ainp,*,iostat=ioerr) porosity,d50,d_90
            if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                          &
               "POROSITY, d50 and D_90",ninp,nout)) return
            else
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               read(ainp,*,iostat=ioerr) Rough,d50
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "ROUGH COEFFICIENT, d50",ninp,nout)) return
               call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
               read(ainp,*,iostat=ioerr) nitersol
               if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,                       &
                  "N° ITERATION SOLID",ninp,nout)) return
         endif
      case default
   endselect
! Assignments 
   if (ncord>0) then
      Med(index)%tipo = tipo
      Med(index)%index = index
      Med(index)%NIterSol = nitersol
      Med(index)%den0 = den0
      Med(index)%eps = eps
      Med(index)%celerita = dsqrt(eps/den0)
      Med(index)%alfaMon = alfaMon
      Med(index)%betaMon = betaMon
      Med(index)%visc = visc
      Med(index)%mumx = viscmx
      Med(index)%cons = cons
      Med(index)%taucri = taucri
      Med(index)%cuin = cuin
      Med(index)%phi = phi
      Med(index)%saturated_medium_flag = saturated_medium_flag
      Med(index)%coes = coes
      Med(index)%limiting_viscosity = limiting_viscosity
      Med(index)%Cs = Cs
      Med(index)%RoughCoef = Rough
      Med(index)%d50 = d50
      Med(index)%gran_vol_frac_max = (1.d0 - porosity)
      Med(index)%d_90 = d_90
      Med(index)%codif = codif
      Med(index)%SettlingCoef = Settling
      Med(index)%gamma = gamma
      Med(index)%InitialIntEn = InitialIntEn
      if (nout>0) then
         write(nout,"(1x,a,i3,1x,a)")  "Medium:.....................",         &
            Med(index)%index,"("//Med(index)%tipo//")"
         write(nout,"(1x,a,1p,e12.4)") "Density:....................",         &
            Med(index)%den0
         write(nout,"(1x,a,1p,e12.4)") "Comprimibility:.............",         &
            Med(index)%eps
         write(nout,"(1x,a,1p,e12.4)") "Celerity:...................",         &
            Med(index)%celerita
         write(nout,"(1x,a,1p,e12.4)") "Alpha Monaghan Coeff.:......",         &
            Med(index)%alfaMon 
         write(nout,"(1x,a,1p,e12.4)") "Beta Monaghan Coeff.:.......",         &
            Med(index)%betaMon
         write(nout,"(1x,a,1p,e12.4)") "Dynamic Viscosity:..........",         &
            Med(index)%visc
         write(nout,"(1x,a,1p,e12.4)") "Max. Dynamic Viscosity:.....",         &
            Med(index)%mumx
         write(nout,"(1x,a,1p,e12.4)") "Tau Critical:...............",         &
            Med(index)%taucri
         write(nout,"(1x,a,1p,e12.4)") "Index curve:................",         &
            Med(index)%cuin
         write(nout,"(1x,a,1p,e12.4)") "Friction Angle:.............",         &
            Med(index)%phi
         write(nout,"(1x,a,1p,l12)")   "Saturated medium flag:......",         &
            Med(index)%saturated_medium_flag
         write(nout,"(1x,a,1p,e12.4)") "Cohesion:...................",         &
            Med(index)%coes
         write(nout,"(1x,a,1p,e12.4)") "Limiting viscosity:.........",         &
            Med(index)%limiting_viscosity
         write(nout,"(1x,a,1p,e12.4)") "Consistency:................",         &
            Med(index)%cons
         write(nout,"(1x,a,1p,e12.4)") "Smagorinsky Constant:.......",         &
            Med(index)%Cs
         write(nout,"(1x,a,1p,e12.4)") "Roughness Coeff.:...........",         &
            Med(index)%RoughCoef
         write(nout,"(1x,a,1p,e12.4)") "d50:........................",         &
            Med(index)%d50
         write(nout,"(1x,a,1p,e12.4)") "d90:........................",         &
            Med(index)%d_90
         write(nout,"(1x,a,1p,e12.4)") "Volume fraction (solid):....",         &
            Med(index)%gran_vol_frac_max
         write(nout,"(1x,a,1p,e12.4)") "Diffusion Coeff.:...........",         &
            Med(index)%codif       
         write(nout,"(1x,a,1p,e12.4)") "Settling Velocity Coeff.:...",         & 
            Med(index)%SettlingCoef       
         write(nout,"(1x,a,1p,e12.4)") "Explosion Gamma Coeff.:.....",         &
            Med(index)%Gamma       
         write(nout,"(1x,a,1p,e12.4)") "Initial Internal Energy.:...",         &
            Med(index)%InitialIntEn       
         write(nout,"(1x,a)")  " "
      endif
      Med(index)%visc = visc / den0
      Med(index)%numx = viscmx / den0
! From degrees to radians 
      Med(index)%phi = Med(index)%phi * PIGRECO / 180.  
   endif
   call ReadRiga(ainp,comment,nrighe,ioerr,ninp)
   if (.NOT.ReadCheck(ioerr,ier,nrighe,ainp,"MEDIUM DATA",ninp,nout)) return
enddo
!------------------------
! Deallocations
!------------------------
return
end subroutine ReadInputMedium

