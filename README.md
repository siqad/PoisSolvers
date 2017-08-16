To be used as a collection of numerical solvers of the Poisson's equation
using various implementations.

PoisFFT
  Uses the fftw library to implement a Fourier transform based solver. This solution is used to verify the simple case for poissonsolver below.

poissonsolver
  Uses the Successive Over-relaxation method to solve a generalised Poisson's equation. The generalised equation takes into account variable permittivity, and allows for electrodes to be arbitrarily placed within the solving range.

# PoisSolvers

A C++ solver for the 3D Poisson's equation using the PoisFFT implementation for the Fourier transform

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

#### PoisFFT

You will need to install PoisFFT for the Fourier transform implementation. Visit their github page and follow the installation instructions at

```
https://github.com/LadaF/PoisFFT
```

### Setting up

1. Install PoisFFT
2. Clone this repository to your preferred directory.

Done!

## Description

The solver an iterative solver that uses the FFT to complete its iterations. The solver allows for electrodes to be defined within the sample body, essentially setting an internal boundary condition at the defined electrode surface. The FFT is able to handle boundary conditions at the sample edges, but cannot easily handle internal boundary conditions. The solver overcomes this by going through an iterative algorithm, and checking the potential at the electrode at each iteration. It then tweaks the volumetric charge density at the electrode surfaces until the potentials at electrode surfaces are within defined error ranges. The solver continues on tweaking the volumetric charge density until it arrives at a sufficiently correct 3D potential.

## Running the tool

In testpoisson.cc, changes can be made to customize the problem to be solved by the solver.

### Setting electrodes

Electrodes can be set by calling

```
Electrodes electrodenamehere(xmin, xmax, ymin, ymax, zmin, zmax, potential, workfunction)
```

Remember to draw the electrode into electrodemap with

```
electrodenamehere.draw(ns, ds, Ls, RHS, electrodemap, chi)
```

### Setting physical parameters

Physical parameters such as permittivity and the volumetric charge density can be set by pointing eps or RHS respectively towards to data set.
Please ensure that the data is of size (ns[0]*ns[1]*ns[2]).
An initial guess for the potential can be set by pointing arr to the intial assumption data set.

###

in the PoisFFT/src/ directory, do:

```
./compileandrun.sh
```

This will compile the code using g++, and run the resulting executable file. 2D plots of the potential and the associated charge density will also be produced using GNUplot.

## Authors

* **Nathan Chiu**

## Acknowledgments

*  Vladimir Fuka: PoisFFT - a free parallel fast Poisson solver, Applied Mathematics and Computation, 2015
available online: http://dx.doi.org/10.1016/j.amc.2015.03.011
preprint: http://arxiv.org/abs/1409.8116
