CXX=g++
CXX_FLAGS = -std=c++11 -Wall -Wextra -g

all: poissolver

poissolver: problem.o phys_engine.o electrodes.o poissolver.o poismain.o
	$(CXX) $(CXXFLAGS) -o bin/gcc/poissolver -fopenmp \
	bin/gcc/precisions.o \
	bin/gcc/parameters.o \
	bin/gcc/customfftw3.o \
	bin/gcc/custompfft.o \
	bin/gcc/fft.o \
	bin/gcc/poisfft.o \
	bin/gcc/main.o \
	bin/gcc/f_mpi_comm_c2f.o \
	bin/gcc/c_binding.o \
	bin/gcc/c_binding_c.o \
	bin/gcc/problem.o \
	bin/gcc/phys_engine.o \
	bin/gcc/electrodes.o \
	bin/gcc/poissolver.o \
	bin/gcc/poismain.o \
	-lm -lfftw3 -lfftw3f -lfftw3_omp -lprofiler -lgfortran

poismain.o: src/poismain.cpp src/poissolver.h src/poisfft.h src/phys_engine.h src/electrodes.h
	$(CXX) $(CXX_FLAGS) -o bin/gcc/poismain.o -c -O3 src/poismain.cpp

poissolver.o: src/poissolver.cpp src/poissolver.h
	$(CXX) $(CXX_FLAGS) -o bin/gcc/poissolver.o -c -O3 src/poissolver.cpp

electrodes.o: src/electrodes.cpp src/electrodes.h
	$(CXX) $(CXX_FLAGS) -o bin/gcc/electrodes.o -c -O3 src/electrodes.cpp

phys_engine.o: src/phys_engine.cpp src/phys_engine.h
	$(CXX) $(CXX_FLAGS) -o bin/gcc/phys_engine.o -c -O3 src/phys_engine.cpp

problem.o: src/problem.cpp src/problem.h
	$(CXX) $(CXX_FLAGS) -o bin/gcc/problem.o -c -O3 src/problem.cpp

# clean:
#     rm *o hello
