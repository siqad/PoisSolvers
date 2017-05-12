#include <iostream>
#include "./src/solver.h"

int main() {
    poisson1DJacobi();
    poisson1DGaussSeidel();
    poisson1DSOR();
    return 0;
}
