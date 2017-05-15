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
    s.set_val(s.h2, Nx*Nx/Lx/Lx, Ny*Ny/Ly/Ly, Nz*Nz/Lz/Lz);
    s.set_val( s.solvemethod, SOR);
    s.init_val( s.V, 0);
    s.init_val( s.rho, 1.69e-19/EPS0);
    s.set_BCs(0, 0, 0, 0, 0, 0);

    //Sanity Check
    std::cout << s.N[0] << " " << s.N[1] << " " << s.N[2] << std::endl;
    std::cout << s.L[0] << " " << s.L[1] << " " << s.L[2] << std::endl;
    std::cout << s.h2[0] << " " << s.h2[1] << " " << s.h2[2] << std::endl;
    std::cout << s.solvemethod << std::endl;
    std::cout << s.V.size() << std::endl;

    //call solver
    s.solve();

    return 0;
}
