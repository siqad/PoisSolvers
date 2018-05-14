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

## Description

The solver uses the FEM to solve for the 3D electric potential given an arrangement of electrodes. The solver is meant to be used in conjunction with another tool, and may take some effort to get working as a stand-alone application.
Boundary conditions at the simulation edges are handled using Robin boundary conditions, while the air-material boundary is implicitly handled by aligning the mesh face to the boundary interface.

## Running the tool

The solver is meant to be called from another tool, with the appropriate input and output arguments.

## Author

**Nathan Chiu**
