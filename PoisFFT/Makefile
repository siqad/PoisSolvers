CXX=g++
CXX_FLAGS = -std=c++11 -Wall -Wextra -g

all: poissolver

poissolver: bin/objs/phys_connector.o bin/objs/problem.o bin/objs/phys_engine.o bin/objs/electrodes.o bin/objs/poissolver.o bin/objs/poismain.o
	$(CXX) $(CXXFLAGS) -o bin/objs/poissolver -fopenmp \
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
	bin/objs/phys_connector.o \
	bin/objs/problem.o \
	bin/objs/phys_engine.o \
	bin/objs/electrodes.o \
	bin/objs/poissolver.o \
	bin/objs/poismain.o \
	-lm -lfftw3 -lfftw3f -lfftw3_omp -lprofiler -lgfortran

bin/objs/poismain.o: src/poismain.cpp src/poissolver.h src/poisfft.h src/phys_engine.h src/electrodes.h
	$(CXX) $(CXX_FLAGS) -o bin/objs/poismain.o -c -O3 src/poismain.cpp

bin/objs/poissolver.o: src/poissolver.cpp src/poissolver.h
	$(CXX) $(CXX_FLAGS) -o bin/objs/poissolver.o -c -O3 src/poissolver.cpp

bin/objs/electrodes.o: src/electrodes.cpp src/electrodes.h
	$(CXX) $(CXX_FLAGS) -o bin/objs/electrodes.o -c -O3 src/electrodes.cpp

bin/objs/phys_engine.o: src/phys_engine.cpp src/phys_engine.h
	$(CXX) $(CXX_FLAGS) -o bin/objs/phys_engine.o -c -O3 src/phys_engine.cpp

bin/objs/problem.o: src/problem.cpp src/problem.h
	$(CXX) $(CXX_FLAGS) -o bin/objs/problem.o -c -O3 src/problem.cpp

bin/objs/phys_connector.o: src/phys_connector.cpp src/phys_connector.h
	mkdir -p bin/objs
	$(CXX) $(CXX_FLAGS) -o bin/objs/phys_connector.o -c -O3 src/phys_connector.cpp

clean:
	rm bin/objs/*
