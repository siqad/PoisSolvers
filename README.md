Holds the 3D FEM Poisson solver.

# PoisSolvers

A Python solver for the 3D Poisson's equation using the FEniCS finite element library.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

#### FEniCS

FEniCS is used to solve the FEM problem posed by the user design.

You can install FEniCS by running:
```
sudo apt install python3-dolfin
```

### Running PoisSolvers on Windows

One of the dependencies used in this project, FEniCS, does not have Windows support. If you have no access to machines running Linux or macOS, the only ways to use FEniCS on Windows is either through [Docker](https://www.docker.com/) or [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/about)(WSL).

If your Windows version supports WSL (Windows 10 version 1607 or above, 64-bit), we recommend using WSL as it is more integrated with Windows 10.

#### Installing WSL and dependent packages

These steps were tested on Windows 10 Education Version 1803 with a Ubuntu 16.04.4 LTS running in WSL.

1. Please follow [this guide](https://docs.microsoft.com/en-us/windows/wsl/install-win10) to install WSL. **Choose Ubuntu as your distribution**.

2. After installing WSL on your system, upgrade packages on the system:
```
sudo apt update && sudo apt upgrade
```

3. Install FEniCS (commands taken from [FEniCS' website](https://fenicsproject.org/download/)):
```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install --no-install-recommends fenics
sudo apt-get dist-upgrade
```

4. Install other dependent packages:
```
sudo apt-get install gmsh python3-tk
```

**Still doesn't work, just putting this documentation up to be completed in the future.**

## Description

The solver uses the FEM to solve for the 3D electric potential given an arrangement of electrodes. The solver is meant to be used in conjunction with another tool, and may take some effort to get working as a stand-alone application.
Boundary conditions at the simulation edges are handled using Robin boundary conditions, while the air-material boundary is implicitly handled by aligning the mesh face to the boundary interface.

## Running the tool

The solver is meant to be called from another tool, with the appropriate input and output arguments.

## Author

**Nathan Chiu**
