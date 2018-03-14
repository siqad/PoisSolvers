//
// Created by nathan on 03/05/17.
//
#ifndef POISSONSOLVER_SOLVER_H
#define POISSONSOLVER_SOLVER_H

#include <vector>
#include "electrodes.h"
#include <string>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#define EVEN 0
#define ODD 1
#define JACOBI 1
#define GAUSS_SEIDEL 2
#define SOR 3
#define SOR_GEN 4
#define SOR_BLAS 5
#define MG 6
#define EPS0 8.85418782e-12
#define Q_E 1.6e-19
#define LATTICECONSTANT 1e-3
#define FILENAMESOR_GEN "outfileSOR_GEN"
#define FILENAMESOR "outfileSOR"
#define FILENAMEJAC "outfileJAC"
#define FILENAMEGAU "outfileGAU"
#define FILENAMERHO "outfileRHO.txt"
#define FILENAMEEPS "outfileEPS.txt"
#define FILENAMEELEC "outfileELEC.txt"
#define FILENAMEBLAS "outfileBLAS"
#define FILENAMEMG "outfileMG"
#define MAXERROR 1e-3
#define PI 3.14159265358979323846
#define DIRICHLET 0
#define NEUMANN 1
#define PERIODIC 2
//from http://hyperphysics.phy-astr.gsu.edu/hbase/Tables/photoelec.html
#define WF_GOLD 5.1 //workfunction for gold in eV
#define WF_COPPER 4.7 //workfunction for copper in eV
#define WF_ZINC 4.3 //workfunction for zinc in eV
#define WF_CESIUM 2.1 //workfunction for cesium in eV
#define WF_NICKEL 5.01
//from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
#define CHI_SI 4.05//electron affinity for silicon in eV

void poisson1DJacobi(void);
void poisson1DGaussSeidel(void);
void poisson1DSOR(void);

class Solver {
public:
    //variables
    Solver();
    ~Solver();
    int*  N; //Lattice points in each dimension. N[0] for x, N[1] for y, N[2] for z
    int solvemethod; //1 for Jacobi, 2 for Gauss-Seidel, 3 for SOR
    int boundarytype; //0 for Dirichlet and 1 for Neumann.
    double* V;
    double* rho; //Volume charge density
    double* L; //Physical lengths for each dimension
    double* eps;
    double h2;
    double h;
    std::pair<double, double>* electrodemap; //electrode mapping - 0 for bulk material, X for electrode at with X voltage.

    //functions
    void relaxmg( int, bool*, bool*, bool*, double*, double**, int, int, double*, double*, int, int, int, double* );
    void calc_Neumann(int, int, int);
    void calc_Periodic(int, int, int, double*);
    void set_N(int, int, int);
    void set_L(double, double, double);
    void set_val(int, int, int, int* );
    double* set_val(double, double);
    void set_val(double, double, double, double*);
    double* init_val( double, double*);
    std::pair<double, double>* init_val( double, std::pair<double, double>*);
    void init_rho(void);
    void initEPS(void);
    void solve(void);
    void poisson3DSOR_BLAS(void);
    void poisson3DSOR_gen(void);
    void poisson3Dmultigrid(void);
    void poisson3DSOR(void);
    void poisson3DJacobi(void);
    void poisson3DGaussSeidel(void);
    void set_BCs(double, double, double, double, double, double);
    void write(std::vector<double>&, std::string );
    void write_2D(std::vector<double>&, std::string);
    void write_2D( double* vals, std::string filename);
    void get_a(double *, double *, int);
    void check_eps(double *, bool*);
    void del( int* );
    void del( double* );
    void del( std::pair<double, double>* );
    void create_a( double** );
    void check_exterior( bool*, int, int, int );
    void check_elec( bool* );
    double ohmic_contact( unsigned long int );
    double normal_eps( unsigned long int, double*, double*);
    double normal( unsigned long int, double*);
    void check_diag( boost::numeric::ublas::compressed_matrix<double> );
    void save_voltage( boost::numeric::ublas::vector<double> );
    void mgrestriction( double*, double* );
    void mgprolongation( double*, double* );
};

#endif //POISSONSOLVER_SOLVER_H
