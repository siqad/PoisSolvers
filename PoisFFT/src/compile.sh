#!/bin/sh -x
#Compiles, runs, and plots the resulting file from the FFT.
g++ -std=c++11 -o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/electrodes.o -c -O3 -g electrodes.cpp
g++ -std=c++11 -o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/poissolver.o -c -O3 -g poissolver.cpp
g++ -std=c++11 -o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/poissolver -fopenmp /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/precisions.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/parameters.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/customfftw3.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/custompfft.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/fft.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/poisfft.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/main.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/f_mpi_comm_c2f.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/c_binding.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/c_binding_c.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/electrodes.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/poissolver.o -lm -lfftw3 -lfftw3f -lfftw3_omp -lprofiler -lgfortran