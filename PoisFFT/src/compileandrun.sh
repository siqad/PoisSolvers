#!/bin/sh -x
#Compiles, runs, and plots the resulting file from the FFT.
g++ -o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/cc_testpoisson.o -c -O3 -g testpoisson.cc
g++ -o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/cc_testpoisson -fopenmp /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/precisions.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/parameters.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/customfftw3.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/custompfft.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/fft.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/poisfft.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/main.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/f_mpi_comm_c2f.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/c_binding.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/c_binding_c.o /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/cc_testpoisson.o -lm -lfftw3 -lfftw3f -lfftw3_omp -lprofiler -lgfortran
/home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/cc_testpoisson
gnuplot /home/nathan/git/PoisSolvers/PoisFFT/bin/gcc/test.gp
