# CompLaB3D

**3D Pore-Scale Reactive Transport Simulator**

CompLaB3D is a Lattice Boltzmann framework for simulating coupled fluid flow, solute transport, and geochemical reactions in 3D pore geometries derived from micro-CT scans or synthetic media.

---

## Features

- **D3Q19 flow solver** - Navier-Stokes via Lattice Boltzmann Method
- **D3Q7 transport solver** - Advection-Diffusion equations
- **Reactive transport** - Operator splitting for geochemical reactions
- **MPI parallelization** - Efficient computation for large 3D domains
- **VTK output** - ParaView-compatible visualization

---

## Repository Structure

```
CompLaB3D-Code/
├── LICENSE
├── README.md
├── CompLaB.xml                # Configuration file
├── index.html                 # Documentation website
├── src/                       # Source code
│   ├── complab3d.cpp
│   ├── complab_functions.hh
│   ├── complab_processors.hh
│   └── defineKinetics.hh
├── build/                     # Build directory (create this)
├── input/                     # Input files (create this)
├── output/                    # Output files (create this)
└── versionControl/            # Palabos library (create this)
```

---

## Installation

### Step 1: Clone the repository

```bash
git clone https://github.com/shahram444/CompLaB3D-Code.git
cd CompLaB3D-Code
```

### Step 2: Create required folders

```bash
mkdir build input output versionControl
```

### Step 3: Install Palabos library

Download Palabos from https://palabos.unige.ch/ and place it in the `versionControl/` folder:

```bash
cd versionControl
# Download and extract Palabos here
```

### Step 4: Build

```bash
cd ../build
cmake ..
make -j4
```

---

## Running Simulation

```bash
cd build
mpirun -np 4 ./complab3d ../CompLaB.xml
```

Output files will be saved in the `output/` folder.

---

## Requirements

- C++11 compiler (GCC 4.8+)
- MPI library (OpenMPI or MPICH)
- CMake 3.10+
- Palabos library
- ParaView (for visualization)

---

## Developers

### 3D Version (CompLaB3D)
**Shahram Asgari**  
shahram.asgari@uga.edu  
Department of Marine Sciences, University of Georgia

**Christof Meile**  
cmeile@uga.edu  
Department of Marine Sciences, University of Georgia

### 2D Version (CompLaB)
**Dr. Heewon Jung**  
hjung@cnu.ac.kr  
Department of Geological Sciences, Chungnam National University

---

## License

GNU Affero General Public License v3.0

Built on [Palabos](https://palabos.unige.ch/) (AGPL v3)
