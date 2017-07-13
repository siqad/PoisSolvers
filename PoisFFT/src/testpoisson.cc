/*
MODIFICATIONS (as required by GPL)
2017 May 11
  -Changes to init_rhs to create non-sinusoidal RHS
  -Added #define FILENAME to define
  -removed call to check_solution from main(). Need exact solution to use check_solution (Need to modify for every test)
  -Added file output to main()
  -Added std::cout statements to check program flow.
  -Added notice to GNU GPL file, as required in section 2c.
  -Changed boundary conditions to DIRICHLET

To build, after installing dependencies (fftw, pfft), go to /src/ and do: sudo scons test
The executable will appear as ../bin/gcc/cc_testpoisson
https://arxiv.org/pdf/1208.0901.pdf
*/

#include <cmath>
#include <iostream>
#include "poisfft.h"
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <utility>
#include <algorithm>

const double pi = 3.14159265358979323846;

#define Q_E 1.6e-19
#define MAXERROR 1e-3
#define IND(i,j,k) (i)*(ns[1]*ns[2])+(j)*(ns[2])+k
#define FILENAME ((char*)"outfile.txt")
#define RHOFILE ((char*) "outrho.txt")

//@@@@@@@@@@@@@@@@ ELECTRODE CLASS DEFINITION

class Electrodes{
  public:
    Electrodes(); //default constructor
    Electrodes(double, double, double, double, double, double, double); //parametrized constructor
    ~Electrodes(); //destructor
    double x[2]; //xmin (x[0]) and xmax (x[1])
    double y[2]; //ymin and ymax
    double z[2]; //zmin and zmax
    double potential;   //pointer after conversion of vector
    void draw(const int[3], const double[3], const double[3], double*, std::pair<bool,double>*);

};

Electrodes::Electrodes( void ){}
Electrodes::Electrodes( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double pot){
  x[0] = xmin;
  x[1] = xmax;
  y[0] = ymin;
  y[1] = ymax;
  z[0] = zmin;
  z[1] = zmax;
  potential = pot;
}

Electrodes::~Electrodes( void ){}

//draw the electrode into an electrode map
void Electrodes::draw(const int ns[3], const double ds[3], const double Ls[3], double* RHS, std::pair<bool,double> *electrodemap){
  int i, j, k;
  std::cout << (int) x[0]/Ls[0]*ns[0] << " " << (int) x[1]/Ls[0]*ns[0] << " " << (int) y[0]/Ls[0]*ns[0] << " " << (int) y[1]/Ls[0]*ns[0]
            << " " << (int) z[0]/Ls[0]*ns[0] << " " << (int) z[1]/Ls[0]*ns[0] << std::endl;
  for(int i = (int) x[0]/Ls[0]*ns[0]; i < (int) x[1]/Ls[0]*ns[0]; i++){ //set RHS 0 inside electrodes
    for(int j = (int) y[0]/Ls[1]*ns[1]; j < (int) y[1]/Ls[1]*ns[1]; j++){
      for(int k = (int) z[0]/Ls[2]*ns[2]; k < (int) z[1]/Ls[2]*ns[2]; k++){
        RHS[i*ns[1]*ns[2] + j*ns[2] + k] = 0;
      }
    }
  }
  for( int iter = 0; iter < 2; iter++){
    i = (int) x[iter]/Ls[0]*ns[0]; //xmin first, then xmax
    for(int j = (int) y[0]/Ls[1]*ns[1]; j <= (int) y[1]/Ls[1]*ns[1]; j++){
      for(int k = (int) z[0]/Ls[2]*ns[2]; k <= (int) z[1]/Ls[2]*ns[2]; k++){
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].first = true; //set true for electrode surface.
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].second = potential; //set electrode potential
      }
    }
  }
  for( int iter = 0; iter < 2; iter++){
    j = (int) y[iter]/Ls[1]*ns[1]; //ymin first, then ymax
    for(int i = (int) x[0]/Ls[0]*ns[0]; i <= (int) x[1]/Ls[0]*ns[0]; i++){
      for(int k = (int) z[0]/Ls[2]*ns[2]; k <= (int) z[1]/Ls[2]*ns[2]; k++){
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].first = true; //set true for electrode surface.
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].second = potential; //set electrode potential
      }
    }
  }
  for( int iter = 0; iter < 2; iter++){
    k = (int) z[iter]/Ls[1]*ns[1]; //zmin
    for(int i = (int) x[0]/Ls[0]*ns[0]; i <= (int) x[1]/Ls[0]*ns[0]; i++){
      for(int j = (int) y[0]/Ls[1]*ns[1]; j <= (int) y[1]/Ls[1]*ns[1]; j++){
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].first = true; //set true for electrode surface.
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].second = potential; //set electrode potential
      }
    }
  }
}

//@@@@@@@@@@@@@@@@@@@@ START OF MAIN

void save_file2D(const int[], const double[], const double[], double*, char[]);
void check_solution(const int[], const double[], const double[], double*);
void init_rhs(const int[], const double[], const double[], double*);
void check_error(const int[], double*, double*, std::pair<bool,double>*);

int main(void){
  std::cout << "Modified by Nathan Chiu. Code is offered as is, with no warranty. See LICENCE.GPL for licence info." << std::endl;
  const double Ls[3] = {10.0, 10.0, 10.0}; //x, y, z domain dimensions
  const int ns[3] = {10, 10, 10}; //x, y, z gridpoint numbers
  double ds[3];  // distances between gridpoints
  const int BCs[6] = {PoisFFT::DIRICHLET, PoisFFT::DIRICHLET,  //boundary conditions
                      PoisFFT::DIRICHLET, PoisFFT::DIRICHLET,
                      PoisFFT::DIRICHLET, PoisFFT::DIRICHLET};
  int i;
  for (i = 0; i<3; i++){ // set the grid, depends on the boundary conditions
    ds[i] = Ls[i] / ns[i];
  }
  double *arr = new double[ns[0]*ns[1]*ns[2]]; // allocate the arrays contiguously, you can use any other class
  double *RHS = new double[ns[0]*ns[1]*ns[2]]; // from which you can get a pointer to contiguous buffer
  double *correction = new double[ns[0]*ns[1]*ns[2]]; // correction used to update RHS
  std::pair<bool,double> *electrodemap = new std::pair<bool,double>[ns[0]*ns[1]*ns[2]]; //stores electrode surface info and potentials.
  init_rhs(ns, ds, Ls, RHS); // set the right-hand side
  Electrodes elec1(3.0, 6.0, 3.0, 6.0, 3.0, 6.0, 10);
  elec1.draw(ns, ds, Ls, RHS, electrodemap);
  PoisFFT::Solver<3, double> S(ns, Ls, BCs); // create solver object, 3 dimensions, double precision
  S.execute(arr, RHS); //run the solver, can be run many times for different right-hand side
  check_error(ns, arr, correction, electrodemap);
  save_file2D(ns, ds, Ls, RHS, RHOFILE);
  save_file2D(ns, ds, Ls, arr, FILENAME); //solution is in arr
  //check_solution(ns, ds, Ls, arr); // check correctness (compares with known, exact solution)
  std::cout << "Ending, deleting variables" << std::endl;
  delete[] RHS;
  delete[] arr;
  delete[] electrodemap;
  delete[] correction;
}

double check_error(const int ns[3], double *arr, double *correction, std::pair<bool,double> *electrodemap){
  double err = 1.1*MAXERROR;
  for(int i = 0; i < ns[0]*ns[1]*ns[2]; i++){
    if(electrodemap[i].first == true){ //only check error at electrodes.
      correction[i] = (electrodemap[i].second-arr[i])/electrodemap[i].second; //intended potential - found potential.
      err = std::max(err, correction[i]); //get largest error value.
    }
  }
  return err;
}


void save_file2D(const int ns[3], const double ds[3], const double Ls[3], double* arr, char fname[]){
  std::ofstream outfile;
  outfile.open(fname, std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping to " << fname << std::endl;
  const int k = ns[2]/2;
  for (int i = 0; i < ns[0]; i++){
    for (int j = 0; j < ns[1]; j++){
        outfile << std::setprecision(5) << std::scientific << i*ds[0] <<
        " " << j*ds[1] << " " << arr[i*ns[1]*ns[2]+j*ns[2]+k] << std::endl;
    }
    outfile << std::endl;
  }
}

void check_solution(const int ns[3], const double ds[3], const double Ls[3], double* a){
  //if the closed form solution is known, then check the obtained solution for correctness
  int i,j,k;
  double sum = 0;
  std::cout << "Checking solution" << std::endl;
  for (i=0;i<ns[0];i++){
    for (j=0;j<ns[1];j++){
      for (k=0;k<ns[2];k++){
        double x = ds[0]*(i+0.5);
        double y = ds[1]*(j+0.5);
        double z = ds[2]*(k+0.5);
        sum += pow(a[IND(i,j,k)] -
            (-1.0/(pow(3*pi/Ls[0],2.0) + pow(5*pi/Ls[1],2.0) + pow(7*pi/Ls[2],2.0)) *
            sin(3*pi*x/Ls[0])*sin(5*pi*y/Ls[1])*sin(7*pi*z/Ls[2])),2.0);
      }
    }
  }
  if (sum < 1e-10) {
    std::cout << "OK, converged\n" << std::endl;
  }
  else {
    std::cout << "FAIL, residuum, did not converge." << sum << std::endl;
  }
}

void init_rhs(const int ns[3], const double ds[3], const double Ls[3], double* a){
  int i,j,k;
  std::cout << "Initialising RHS" << std::endl;
  for (i=0;i<ns[0];i++){
    double x = ds[0]*(i+0.5);
    for (j=0;j<ns[1];j++){
      double y = ds[1]*(j+0.5);
      for (k=0;k<ns[2];k++){
        double z = ds[2]*(k+0.5);
        a[IND(i,j,k)] = 1.5e10*Q_E; //Set default set to bulk volume charge density.
        // setting to 0 inside electrodes dealt with in Electrodes::draw()
        // if( (x >= Ls[0]/3.0 && x <= Ls[0]*2.0/3.0) && (y >= Ls[1]/3.0 && y <= Ls[1]*2.0/3.0) && (z >= Ls[2]/3.0 && z <= Ls[2]*2.0/3.0) ){
        //   a[IND(i,j,k)] = 0.0; //set to 0 inside of electrodes.
        // }
      }
    }
  }
  std::cout << "Finished initialisation" << std::endl;
}
