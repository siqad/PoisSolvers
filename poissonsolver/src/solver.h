//
// Created by nathan on 03/05/17.
//
#ifndef POISSONSOLVER_SOLVER_H
#define POISSONSOLVER_SOLVER_H

#include <vector>

#define JACOBI 1
#define GAUSS_SEIDEL 2
#define SOR 3
#define SOR_GEN 4
#define EPS0 8.85418782e-12
#define Q_E 1.6e-19
#define LATTICECONSTANT 1e-3
#define FILENAMESOR_GEN "outfileSOR_GEN.txt"
#define FILENAMESOR "outfileSOR.txt"
#define FILENAMEJAC "outfileJAC.txt"
#define FILENAMEGAU "outfileGAU.txt"
#define FILENAMERHO "outfileRHO.txt"
#define FILENAMEEPS "outfileEPS.txt"
#define MAXERROR 1e-5
#define IND(i,j,k) (i)*(N[1]*N[2])+(j)*(N[2])+k

void poisson1DJacobi(void);
void poisson1DGaussSeidel(void);
void poisson1DSOR(void);

class Solver {
  public:
    //variables
    int solvemethod; //1 for Jacobi, 2 for Gauss-Seidel, 3 for SOR
    double h2;
    std::vector<int>  N; //Lattice points in each dimension. N[0] for x, N[1] for y, N[2] for z
    std::vector<double>  L; //Physical lengths for each dimension
    std::vector<double>  V; //Potential vector.
//    std::vector<double>  h2; //Lattice spacing (lattice point per length) squared
    std::vector<double>  rho; //Volume charge density
    std::vector<double>  eps; //Relative permittivity
    std::vector<int> temp; //temp variable for developing
    //functions
    void set_N( int, int, int );
    void set_L( double, double, double );
    void set_val( int&, int);
    void set_val( double&, double);
    void set_val( std::vector<double>&, double, double, double);
    void set_val( std::vector<int>&, int, int, int);
    void init_val( std::vector<double>&, double );
    void init_rho( void );
    void init_eps( void );
    void solve( void );
    void poisson3DSOR_gen( void );
    void poisson3DSOR( void );
    void poisson3DJacobi( void );
    void poisson3DGaussSeidel( void );
    void set_BCs( double, double, double, double, double, double);
    void write( std::vector<double>&, std::string );
    std::vector<double> get_a( std::vector<double>&, int);
};

#endif //POISSONSOLVER_SOLVER_H
