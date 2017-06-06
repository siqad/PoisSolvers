To be used as a collection of numerical solvers of the Poisson's equation
using various implementations.

PoisFFT
  Uses the fftw library to implement a Fourier transform based solver. This solution is used to verify the simple case for poissonsolver below.

poissonsolver
  Uses the Successive Over-relaxation method to solve a generalised Poisson's equation. The generalised equation takes into account variable permittivity, and allows for electrodes to be arbitrarily placed within the solving range.
