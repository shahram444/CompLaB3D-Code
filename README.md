# CompLaB3D

**Pore-Scale Reactive Transport Simulator**

CompLaB3D is a Lattice Boltzmann framework for simulating coupled 
fluid flow, solute transport, and geochemical reactions in 3D pore 
geometries derived from micro-CT scans or synthetic media.

🌐 **[Documentation & Website](https://YOUR-USERNAME.github.io/CompLaB3D/)**

## Features

- D3Q19 flow solver (Navier-Stokes)
- D3Q7 transport solver (Advection-Diffusion)
- Reactive transport with operator splitting
- MPI parallelization for large domains
- VTK output for ParaView visualization

## Quick Start
```bash
# Build
mkdir build && cd build
cmake ..
make -j4

# Run
mpirun -np 4 ./complab3d config.xml
```

## Requirements

- C++11 compiler (GCC 4.8+)
- MPI library (OpenMPI or MPICH)
- CMake 3.10+
- ParaView (for visualization)

## Built On

[Palabos](https://palabos.unige.ch/) - Open-source Lattice Boltzmann library

## Acknowledgements

| Role | Name |
|------|------|
| 2D Version | Dr. Heewon Jung |
| 3D Extension | Shahram Asgari |
| PI | Dr. Christof Meile |
| Affiliation | University of Georgia, Dept. of Marine Sciences |

## License

GNU Affero General Public License v3.0

This software is built on Palabos (AGPL v3).
