#include <iostream>
#include "./src/solver.h"

int main() {
/*
    poisson1DJacobi();
    poisson1DGaussSeidel();
    poisson1DSOR();
*/
    Solver s;
    s.set_N(1, 2, 3);
    s.set_L(4, 5, 6);
    s.set_method(SOR);
    s.init_V();
    std::cout << s.N[0] << " " << s.N[1] << " " << s.N[2] << std::endl;
    std::cout << s.L[0] << " " << s.L[1] << " " << s.L[2] << std::endl;
    std::cout << s.solvemethod << std::endl;
    std::cout << s.V.size() << std::endl;

    return 0;
}
