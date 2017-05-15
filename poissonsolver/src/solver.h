//
// Created by nathan on 03/05/17.
//
#ifndef POISSONSOLVER_SOLVER_H
#define POISSONSOLVER_SOLVER_H

#include <vector>

#define JACOBI 1
#define GAUSS_SEIDEL 2
#define SOR 3

void poisson1DJacobi(void);
void poisson1DGaussSeidel(void);
void poisson1DSOR(void);

class Solver {
  public:
    //variables
    int solvemethod; //1 for Jacobi, 2 for Gauss-Seidel, 3 for SOR
    std::vector<int>  N; //Lattice points in each dimension. N[0] for x, N[1] for y, N[2] for z
    std::vector<double>  L; //Physical lengths for each dimension
    std::vector<double>  V; //Potential vector.
    std::vector<double>  h2; //Lattice spacing (length per lattice point) squared
    std::vector<double>  rho; //Volume charge density
    std::vector<int> temp; //temp variable for developing
    //functions
    void set_N( int, int, int );
    void set_L( double, double, double );
    void set_val( int&, int);
    void set_val( std::vector<double>&, double, double, double);
    void set_val( std::vector<int>&, int, int, int);
    void init_val( std::vector<double>&, double );
    void solve( void );
    void poisson3DSOR( std::vector<double> );
};

#endif //POISSONSOLVER_SOLVER_H
