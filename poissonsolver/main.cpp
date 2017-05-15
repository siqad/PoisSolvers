#include <iostream>
#include "./src/solver.h"
int main() {
/*
    poisson1DJacobi();
    poisson1DGaussSeidel();
    poisson1DSOR();
*/

    //Parameters
    int Nx = 1;
    int Ny = 2;
    int Nz = 3;
    double Lx = 4;
    double Ly = 5;
    double Lz = 6;


    //initialize
    Solver s;
    s.set_val(s.N, Nx, Ny, Nz);
    s.set_val(s.L, Lx, Ly, Lz);
    s.set_val(s.h2, Lx*Lx/Nx/Nx, Ly*Ly/Ny/Ny, Lz*Lz/Nz/Nz);
    s.set_val( s.solvemethod, SOR);
    s.init_val( s.V, 0);
    s.init_val( s.rho, 1.69e-19);

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
