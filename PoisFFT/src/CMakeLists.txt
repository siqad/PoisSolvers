cmake_minimum_required(VERSION 3.7)
project(src)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        c_binding.c
        f_mpi_comm_c2f.c
        missing_pfft.c
        poisfft.h
        electrodes.h
        electrodes.cpp
        testpoisson.c
        testpoisson.cc)

add_executable(src ${SOURCE_FILES})
