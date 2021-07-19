# SPHERA
SPHERA v.9.0.0(RSE SpA) is freeresearch software(FOSS) based on the SPH (“Smoothed Particle Hydrodynamics”) method, which represents a mesh-less Computational Fluid Dynamics techniquefor free surface and multi-phase flows. So far, SPHERA has been applied to represent several types of floods (with transport of solid bodies and bed-load transportflood-control works, flood-induced damage; domain spatial coverage of some hundredths of squared kilometres) and fast landslides, sloshing tanks, sea waves and sediment removal from water reservoirs.

## Getting Start
```bash
git clone https://github.com/AndreaAmicarelliRSE/SPHERA.git
cd SPHERA
```
### Build with [fortran-lang/fpm](https://github.com/fortran-lang/fpm)
Fortran Package Manager (fpm) is a great package manager and build system for Fortran.
You can build using provided `fpm.toml`:
```bash
fpm build --flag "-cpp -DSPACE_3D -fallow-invalid-boz"  # user defined flag
```

### Build with make
You can build using provided `Makefile`:
```bash
cd src
make COMPILATION_FLAGS="-fallow-invalid-boz"
```