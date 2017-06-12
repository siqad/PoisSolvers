#ifndef ELECTRODES_H
#define ELECTRODES_H

#include <vector>
#include <string>
#include "solver.h"

class Solver;

class Electrodes{
  public:
    Electrodes();
    ~Electrodes();
    int* pcentre; //pointer after conversion of vector
    int* pdims;   //pointer after conversion of vector
    std::vector<int> centre; //x, y, z for centre of rectangle (lattice space)
    std::vector<int> dims; //x, y, z dimensions for rectangle (in lattice units)
    double potential; //voltage of electrode
    double workfunction; //workfunction of material.
    void init_elec( void );
    void draw( Solver * );
    void convert( void );
};

#endif //ELECTRODES_H
