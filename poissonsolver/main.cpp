#include <iostream>
#include "./src/solver.h"

int main() {
/*
    poisson1DJacobi();
    poisson1DGaussSeidel();
    poisson1DSOR();
*/
    //Parameters
    int Nx = 100;
    int Ny = 100;
    int Nz = 100;
    double Lx = 1;
    double Ly = 2;
    double Lz = 3;
    //initialize
    Solver s;
    s.set_val(s.N, Nx, Ny, Nz);
    s.set_val(s.L, Lx, Ly, Lz);
    s.set_val(s.h2, LATTICECONSTANT*LATTICECONSTANT);
    s.init_val( s.V, 0);
    s.init_val( s.rho, 1.69e-19/EPS0);
    s.set_BCs(0, 0, 0, 0, 0, 0);
    s.set_val( s.solvemethod, SOR);
    //call solver
    s.solve();
    s.write();
    std::cout << s.V[643178] << std::endl;
    return 0;
}
