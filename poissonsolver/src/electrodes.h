#ifndef ELECTRODES_H
#define ELECTRODES_H

#include <vector>
#include <string>
#include "solver.h"

class Solver;

class Electrodes{
  public:
    std::vector<int> centre; //x, y, z for centre of rectangle (lattice space)
    std::vector<int> dims; //x, y, z dimensions for rectangle (in lattice units)
    double workfunc;
    void init_elec( Solver * );
    void draw( Solver * );
};

#endif //ELECTRODES_H
