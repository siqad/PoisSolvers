cmake_minimum_required(VERSION 3.7)
project(poissonsolver_ref)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/cptstd.h
        src/matrix.h
        src/solvers.cpp
        main.cpp)

add_executable(poissonsolver_ref ${SOURCE_FILES})