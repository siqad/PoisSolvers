cmake_minimum_required(VERSION 3.7)
project(poissonsolver)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        src/solver.cpp
        src/solver.h
        src/electrodes.cpp
        src/electrodes.h
        main.cpp)

add_executable(poissonsolver ${SOURCE_FILES})
