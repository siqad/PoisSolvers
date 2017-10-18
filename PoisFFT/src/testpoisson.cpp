/*
MODIFICATIONS
  -Changes to init_rhs to create non-sinusoidal RHS
  -Removed call to check_solution from main(). Need exact solution to use check_solution (Need to modify for every test)
  -Added file output
  -Added Electrode class, constructor, and drawing functionality
  -Added electrodemap, a centralized place to save electrode information
  -Implement algorithm for internal boundary condition (at electrodes)
  -Add workfunction and electron affinity functionality.
  -Add iterative looping, now performs FFT many times, each iteration approaching the desired answer
  -Adjusted weight until iterative loop is stable for variable input parameters.
  -Added google-perftools profiler
  -Added consideration for relative permittivity
Conceptually, FFT gives a solution for given charge density and permittivity. Since it is not clear how to set "internal boundary conditions"
to fix the electrode potentials in the FFT solution, the program uses the FFT to iteratively produce solutions to Poisson's equation,
varying charge at the electrode surfaces until the potential at all electrode surfaces is sufficiently accurate.
*/

#include <cmath>
#include <iostream>
#include "poisfft.h"
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <utility>
#include <algorithm>
#include <gperftools/profiler.h>
#include <ctime>
#include "electrodes.hpp"

const double pi = 3.14159265358979323846;

#define EPS0 8.85418782e-12
#define Q_E 1.6e-19
#define MAXERROR 5e-2
#define IND(i,j,k) (i)*(ns[1]*ns[2])+(j)*(ns[2])+k
#define FILENAME ((char*)"outfile.txt")
#define RHOFILE ((char*) "outrho.txt")
#define EPSFILE ((char*) "outeps.txt")
#define CORRECTIONFILE ((char*) "outcorr.txt")
#define WF_GOLD 5.1 //workfunction for gold in eV
#define WF_COPPER 4.7 //workfunction for copper in eV
#define WF_ZINC 4.3 //workfunction for zinc in eV
#define WF_CESIUM 2.1 //workfunction for cesium in eV
#define WF_NICKEL 5.01
#define CHI_SI 4.05//electron affinity for silicon in eV from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
#define EPS_SI 11.7

// class Electrodes{
//   public:
//     Electrodes(); //default constructor
//     Electrodes(double, double, double, double, double, double, double, double); //parametrized constructor
//     ~Electrodes(); //destructor
//     double x[2]; //xmin (x[0]) and xmax (x[1])
//     double y[2]; //ymin and ymax
//     double z[2]; //zmin and zmax
//     double potential;   //pointer after conversion of vector
//     double WF;
//     void draw(const int[3], const double[3], const double[3], double*, std::pair<int,double>*, double*);
// };
//
// Electrodes::Electrodes( void ){}
// Electrodes::Electrodes( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double pot, double wf){
//   x[0] = xmin;
//   x[1] = xmax;
//   y[0] = ymin;
//   y[1] = ymax;
//   z[0] = zmin;
//   z[1] = zmax;
//   potential = pot;
//   WF = wf;
// }
//
// Electrodes::~Electrodes( void ){}
// //https://ecee.colorado.edu/~bart/book/book/chapter3/ch3_2.htm#3_2_3 for information on metal-semiconductor junction under bias
// void Electrodes::draw(const int ns[3], const double ds[3], const double Ls[3], double* RHS, std::pair<int,double> *electrodemap, double* chi){
//   int i, j, k; //draw the electrode into an electrode map
//   for(i = (int) ns[0]*x[0]/Ls[0]; i < (int) ns[0]*x[1]/Ls[0]; i++){ //set RHS 0 inside electrodes
//     for(j = (int) ns[1]*y[0]/Ls[1]; j < (int) ns[1]*y[1]/Ls[1]; j++){
//       for(k = (int) ns[2]*z[0]/Ls[2]; k < (int) ns[2]*z[1]/Ls[2]; k++){
//         RHS[IND(i,j,k)] = 0;
//       }
//     }
//   }
//   for( int iter = 0; iter < 2; iter++){
//     i = (int) ns[0]*x[iter]/Ls[0]; //xmin first, then xmax
//     for(j = (int) ns[1]*y[0]/Ls[1]; j <= (int) ns[1]*y[1]/Ls[1]; j++){
//       for(k = (int) ns[2]*z[0]/Ls[2]; k <= (int) ns[2]*z[1]/Ls[2]; k++){
//         electrodemap[IND(i,j,k)].first = true; //set true for electrode surface.
//         electrodemap[IND(i,j,k)].second = potential; //set electrode potential
//         if((x[0]!=0 && iter==0) || (x[1]!=ns[0]-1 && iter==1)){ //set potential of adjacent silicon with WF.
//           electrodemap[IND(i-1+2*iter,j,k)].first = true;
//           electrodemap[IND(i-1+2*iter,j,k)].second = potential - (WF-chi[IND(i-1+2*iter,j,k)]);
//         }
//       }
//     }
//     j = (int) ns[1]*y[iter]/Ls[1]; //ymin first, then ymax
//     for(i = (int) ns[0]*x[0]/Ls[0]; i <= (int) ns[0]*x[1]/Ls[0]; i++){
//       for(k = (int) ns[2]*z[0]/Ls[2]; k <= (int) ns[2]*z[1]/Ls[2]; k++){
//         electrodemap[IND(i,j,k)].first = true; //set true for electrode surface.
//         electrodemap[IND(i,j,k)].second = potential; //set electrode potential
//         if((y[0]!=0 && iter==0) || (y[1]!=ns[1]-1 && iter==1)){ //set potential of adjacent silicon with WF.
//           electrodemap[IND(i,j-1+2*iter,k)].first = true;
//           electrodemap[IND(i,j-1+2*iter,k)].second = potential - (WF-chi[IND(i,j-1+2*iter,k)]);
//         }
//       }
//     }
//     k = (int) ns[2]*z[iter]/Ls[2]; //zmin
//     for(i = (int) ns[0]*x[0]/Ls[0]; i <= (int) ns[0]*x[1]/Ls[0]; i++){
//       for(j = (int) ns[1]*y[0]/Ls[1]; j <= (int) ns[1]*y[1]/Ls[1]; j++){
//         electrodemap[IND(i,j,k)].first = true; //set true for electrode surface.
//         electrodemap[IND(i,j,k)].second = potential; //set electrode potential
//         if((z[0]!=0 && iter==0) || (z[1]!=ns[2]-1 && iter==1)){ //set potential of adjacent silicon with WF.
//           electrodemap[IND(i,j,k-1+2*iter)].first = true;
//           electrodemap[IND(i,j,k-1+2*iter)].second = potential - (WF-chi[IND(i,j,k-1+2*iter)]);
//         }
//       }
//     }
//   }
// }

void save_file2D(const int[], const double[], const double[], double*, char[]);
void check_solution(const int[], const double[], const double[], double*);
void init_eps(const int[], const double[], const double[], double*);
void init_rhs(const int[], const double[], const double[], double*, double*, double*);
double check_error(const int[], const double[], double*, double*, std::pair<int,double>*, int*, double*, double*, const int[]);
void apply_correction(const int[], const double[], double*, double*, std::pair<int,double>*);

int main(void){
  std::cout << "Modified by Nathan Chiu. Code is offered as is, with no warranty. See LICENCE.GPL for licence info." << std::endl;
  // const double Ls[3] = {10.0, 10.0, 10.0}; //x, y, z domain dimensions in DECAMETRE
  // const double Ls[3] = {1.0, 1.0, 1.0}; //x, y, z domain dimensions in METRE
  // const double Ls[3] = {1.0e-1, 1.0e-1, 1.0e-1}; //x, y, z domain dimensions in DECIMETRE
  const double Ls[3] = {1.0e-6, 1.0e-6, 1.0e-6}; //x, y, z domain dimensions in CENTIMETRE
  // const double Ls[3] = {1.0e-3, 1.0e-3, 1.0e-3}; //x, y, z domain dimensions in MILLIMETRE
  // const double Ls[3] = {1.0e-6, 1.0e-6, 1.0e-6}; //x, y, z domain dimensions in MICROMETRE
  // const double Ls[3] = {1.0e-9, 1.0e-9, 1.0e-9}; //x, y, z domain dimensions in NANOMETRE
  const int ns[3] = {150, 150, 150}; //x, y, z gridpoint numbers
  double ds[3];  // distances between gridpoints
  double cycleErr;
  int* indexErr = new int;
  const int BCs[6] = {PoisFFT::PERIODIC, PoisFFT::PERIODIC,  //boundary conditions
                      PoisFFT::PERIODIC, PoisFFT::PERIODIC,
                      PoisFFT::PERIODIC, PoisFFT::PERIODIC};
  int i;
  for (i = 0; i<3; i++){ // set the grid, depends on the boundary conditions
    ds[i] = Ls[i] / ns[i];
  }
  double *arr = new double[ns[0]*ns[1]*ns[2]]; // allocate the arrays contiguously, you can use any other class
  double *arrOld = new double[ns[0]*ns[1]*ns[2]]; // allocate the arrays contiguously, you can use any other class
  double *eps = new double[ns[0]*ns[1]*ns[2]]; //RELATIVE permittivity.
  double *chi = new double[ns[0]*ns[1]*ns[2]]; //electron affinity or WF
  double *RHS = new double[ns[0]*ns[1]*ns[2]]; // from which you can get a pointer to contiguous buffer, will contain rho.
  double *correction = new double[ns[0]*ns[1]*ns[2]]; // correction used to update RHS
  std::pair<int,double> *electrodemap = new std::pair<int,double>[ns[0]*ns[1]*ns[2]]; //stores electrode surface info and potentials.
  int cycleCount = 0;
  const std::clock_t begin_time = std::clock();
  init_eps(ns, ds, Ls, eps); // set permittivity
  init_rhs(ns, ds, Ls, chi, eps, RHS); // set rho and apply relative permittivity, also save electron affinity for bulk.
  Electrodes elec1(0.2e-6, 0.4e-6, 0.2e-6, 0.4e-6, 0.3e-6, 0.6e-6, 0, WF_GOLD);
  Electrodes elec2(0.2e-6, 0.4e-6, 0.6e-6, 0.8e-6, 0.3e-6, 0.6e-6, 5, WF_GOLD);
  Electrodes elec3(0.6e-6, 0.8e-6, 0.2e-6, 0.4e-6, 0.3e-6, 0.6e-6, 10, WF_GOLD);
  Electrodes elec4(0.6e-6, 0.8e-6, 0.6e-6, 0.8e-6, 0.3e-6, 0.6e-6, 15, WF_GOLD);
  elec1.draw(ns, ds, Ls, RHS, electrodemap, chi); //separately call draw for each electrode.
  elec2.draw(ns, ds, Ls, RHS, electrodemap, chi);
  elec3.draw(ns, ds, Ls, RHS, electrodemap, chi);
  elec4.draw(ns, ds, Ls, RHS, electrodemap, chi);
  PoisFFT::Solver<3, double> S(ns, Ls, BCs); //   create solver object, 3 dimensions, double precision
  std::cout << "Beginning solver" << std::endl;
  ProfilerStart("./profileresult.out"); //using google performance tools
  do{
    S.execute(arr, RHS); //run the solver, can be run many times for different right-hand side
    cycleErr = check_error(ns, Ls, arr, correction, electrodemap, indexErr, arrOld, eps, BCs);
    apply_correction(ns, ds, RHS, correction, electrodemap);
    std::cout << "On cycle " << cycleCount << " with error " << cycleErr << " at index " << *indexErr << ". " << arr[*indexErr] << " " << electrodemap[*indexErr].second << std::endl;
    cycleCount++;
  }while(cycleErr > MAXERROR && cycleErr != 0);
  ProfilerStop();
  std::cout << "Finished on cycle " << cycleCount << " with error " << cycleErr << std::endl;
  save_file2D(ns, ds, Ls, RHS, RHOFILE);
  save_file2D(ns, ds, Ls, arr, FILENAME); //solution is in arr
  std::cout << "Time elapsed: " << float(clock()-begin_time)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "Ending, deleting variables" << std::endl;
  delete indexErr;
  delete[] eps;
  delete[] RHS;
  delete[] arr;
  delete[] electrodemap;
  delete[] chi;
  delete[] correction;
}

void apply_correction(const int ns[3], const double ds[3], double *RHS, double *correction, std::pair<int,double> *electrodemap){
  for(int i = 0; i < ns[0]*ns[1]*ns[2]; i++){
    if(electrodemap[i].first == true){ //only correct error at electrode surfaces.
      RHS[i] -= correction[i];
    }
  }
}

double check_error(const int ns[3], const double Ls[3], double *arr, double *correction, std::pair<int,double> *electrodemap, int *indexErr, double *arrOld, double *eps, const int BCs[6]){
  double err = 0;
  double errOld;
  double correctionWeight;
  if(BCs[0] == PoisFFT::PERIODIC){
    correctionWeight = 1e-7*ns[0]*EPS0/Q_E/Ls[0]/Ls[0]; //Periodic
  } else if(BCs[0] == PoisFFT::DIRICHLET || BCs[0] == PoisFFT::NEUMANN){
    correctionWeight = 0.5e-7*ns[0]*EPS0/Q_E/Ls[0]/Ls[0];  //Dirichlet & Neumann
  }
  for(int i = 0; i < ns[0]*ns[1]*ns[2]; i++){
    if(electrodemap[i].first == true){ //only check error at electrodes.
      errOld = err;
      correction[i] = electrodemap[i].second-arr[i]; //intended potential - found potential. Looking at laplacian of potential doesn't allow snapping to electrode potentials.
      correction[i] *= correctionWeight;
      if(electrodemap[i].second != 0){
        err = std::max(err, fabs((arr[i] - electrodemap[i].second)/electrodemap[i].second)); //get largest error value.
      } else {
        err = std::max(err, arr[i]);
      }
      if( errOld != err ){
        *indexErr = i;
      }
    }
    arrOld[i] = arr[i];
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
        outfile << std::setprecision(5) << std::scientific << i*ds[0] << " " << j*ds[1] << " " << arr[IND(i,j,k)] << std::endl;
    }
    outfile << std::endl;
  }
}

void init_eps(const int ns[3], const double ds[3], const double Ls[3], double* a){
  int i,j,k;
  std::cout << "Initialising rho" << std::endl;
  for (i=0;i<ns[0];i++){
    double x = ds[0]*(i+0.5);
    for (j=0;j<ns[1];j++){
      double y = ds[1]*(j+0.5);
      for (k=0;k<ns[2];k++){
        double z = ds[2]*(k+0.5);
        if (y < Ls[1]/2){
          a[IND(i,j,k)] = EPS_SI; //Si relative permittivity
          // a[IND(i,j,k)] = 1;  //Free space
        } else{
          // a[IND(i,j,k)] = EPS_SI; //Si relative permittivity
          a[IND(i,j,k)] = 1;  //Free space
        }
      }
    }
  }
  std::cout << "Finished rho initialisation" << std::endl;
}

void init_rhs(const int ns[3], const double ds[3], const double Ls[3], double* chi, double* eps, double* a){
  int i,j,k;
  std::cout << "Initialising RHS" << std::endl;
  for (i=0;i<ns[0];i++){
    double x = ds[0]*(i+0.5);
    for (j=0;j<ns[1];j++){
      double y = ds[1]*(j+0.5);
      for (k=0;k<ns[2];k++){
        double z = ds[2]*(k+0.5);
        a[IND(i,j,k)] = 1e16*Q_E/EPS0/eps[IND(i,j,k)]; //in m^-3, scale by permittivity
        chi[IND(i,j,k)] = CHI_SI;
      }
    }
  }
  std::cout << "Finished RHS initialisation" << std::endl;
}
