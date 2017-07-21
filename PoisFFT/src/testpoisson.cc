/*
MODIFICATIONS (as required by GPL)
2017 May 11
  -Changes to init_rhs to create non-sinusoidal RHS
  -Added #define FILENAME to define
  -removed call to check_solution from main(). Need exact solution to use check_solution (Need to modify for every test)
  -Added file output to main()
  -Added std::cout statements to check program flow.
  -Added notice to GNU GPL file, as required in section 2c.
  -Added Electrode class and drawing functionality
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
#include <gperftools/profiler.h>
#include <ctime>

const double pi = 3.14159265358979323846;

#define EPS0 8.85418782e-12
#define Q_E 1.6e-19
#define MAXERROR 1e-3
#define IND(i,j,k) (i)*(ns[1]*ns[2])+(j)*(ns[2])+k
#define FILENAME ((char*)"outfile.txt")
#define RHOFILE ((char*) "outrho.txt")
#define CORRECTIONFILE ((char*) "outcorr.txt")
#define WF_GOLD 5.1 //workfunction for gold in eV
#define WF_COPPER 4.7 //workfunction for copper in eV
#define WF_ZINC 4.3 //workfunction for zinc in eV
#define WF_CESIUM 2.1 //workfunction for cesium in eV
#define WF_NICKEL 5.01
#define CHI_SI 4.05//electron affinity for silicon in eV from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html



class Electrodes{
  public:
    Electrodes(); //default constructor
    Electrodes(double, double, double, double, double, double, double, double); //parametrized constructor
    ~Electrodes(); //destructor
    double x[2]; //xmin (x[0]) and xmax (x[1])
    double y[2]; //ymin and ymax
    double z[2]; //zmin and zmax
    double potential;   //pointer after conversion of vector
    double WF;
    void draw(const int[3], const double[3], const double[3], double*, std::pair<int,double>*, double*);
};

Electrodes::Electrodes( void ){}
Electrodes::Electrodes( double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double pot, double wf){
  x[0] = xmin;
  x[1] = xmax;
  y[0] = ymin;
  y[1] = ymax;
  z[0] = zmin;
  z[1] = zmax;
  potential = pot;
  WF = wf;
}

Electrodes::~Electrodes( void ){}

void Electrodes::draw(const int ns[3], const double ds[3], const double Ls[3], double* RHS, std::pair<int,double> *electrodemap, double* chi){
  int i, j, k; //draw the electrode into an electrode map
  //http://in.ncu.edu.tw/ncume_ee/SchottkyDiode.htm says that potential across metal-semi conductor
  //follows Schottky Barrier model, regardless of forward/reverse bias.
  for(i = (int) ns[0]*x[0]/Ls[0]; i < (int) ns[0]*x[1]/Ls[0]; i++){ //set RHS 0 inside electrodes
    for(j = (int) ns[1]*y[0]/Ls[1]; j < (int) ns[1]*y[1]/Ls[1]; j++){
      for(k = (int) ns[2]*z[0]/Ls[2]; k < (int) ns[2]*z[1]/Ls[2]; k++){
        RHS[IND(i,j,k)] = 0;
        chi[IND(i,j,k)] = WF;
      }
    }
  }
  std::cout << ns[0]*x[1]/Ls[0] << " " << ns[0]-1 << std::endl;
  for( int iter = 0; iter < 2; iter++){
    i = (int) ns[0]*x[iter]/Ls[0]; //xmin first, then xmax
    for(j = (int) ns[1]*y[0]/Ls[1]; j <= (int) ns[1]*y[1]/Ls[1]; j++){
      for(k = (int) ns[2]*z[0]/Ls[2]; k <= (int) ns[2]*z[1]/Ls[2]; k++){
        electrodemap[IND(i,j,k)].first = true; //set true for electrode surface.
        electrodemap[IND(i,j,k)].second = potential; //set electrode potential
        if((x[0]!=0 && iter==0) || (x[1]!=ns[0]-1 && iter==1)){ //set potential of adjacent silicon with WF.
          electrodemap[IND(i-1+2*iter,j,k)].first = true;
          electrodemap[IND(i-1+2*iter,j,k)].second = potential - (WF-chi[IND(i-1+2*iter,j,k)]);
        }
      }
    }
    j = (int) ns[1]*y[iter]/Ls[1]; //ymin first, then ymax
    for(i = (int) ns[0]*x[0]/Ls[0]; i <= (int) ns[0]*x[1]/Ls[0]; i++){
      for(k = (int) ns[2]*z[0]/Ls[2]; k <= (int) ns[2]*z[1]/Ls[2]; k++){
        electrodemap[IND(i,j,k)].first = true; //set true for electrode surface.
        electrodemap[IND(i,j,k)].second = potential; //set electrode potential
        if((y[0]!=0 && iter==0) || (y[1]!=ns[1]-1 && iter==1)){ //set potential of adjacent silicon with WF.
          electrodemap[IND(i,j-1+2*iter,k)].first = true;
          electrodemap[IND(i,j-1+2*iter,k)].second = potential - (WF-chi[IND(i,j-1+2*iter,k)]);
        }
      }
    }
    k = (int) ns[2]*z[iter]/Ls[2]; //zmin
    for(i = (int) ns[0]*x[0]/Ls[0]; i <= (int) ns[0]*x[1]/Ls[0]; i++){
      for(j = (int) ns[1]*y[0]/Ls[1]; j <= (int) ns[1]*y[1]/Ls[1]; j++){
        electrodemap[IND(i,j,k)].first = true; //set true for electrode surface.
        electrodemap[IND(i,j,k)].second = potential; //set electrode potential
        if((z[0]!=0 && iter==0) || (z[1]!=ns[2]-1 && iter==1)){ //set potential of adjacent silicon with WF.
          electrodemap[IND(i,j,k-1+2*iter)].first = true;
          electrodemap[IND(i,j,k-1+2*iter)].second = potential - (WF-chi[IND(i,j,k-1+2*iter)]);
        }
      }
    }
  }
}

void save_file2D(const int[], const double[], const double[], double*, char[]);
void check_solution(const int[], const double[], const double[], double*);
void init_eps(const int[], const double[], const double[], double*);
void init_rhs(const int[], const double[], const double[], double*, double*, double*);
double check_error(const int[], const double[], double*, double*, std::pair<int,double>*, double*, const int[]);
void apply_correction(const int[], const double[], double*, double*, std::pair<int,double>*);

int main(void){
  std::cout << "Modified by Nathan Chiu. Code is offered as is, with no warranty. See LICENCE.GPL for licence info." << std::endl;
  // const double Ls[3] = {0.1, 0.1, 0.1}; //x, y, z domain dimensions
  const double Ls[3] = {1.0, 1.0, 1.0}; //x, y, z domain dimensions
  // const double Ls[3] = {10.0, 10.0, 10.0}; //x, y, z domain dimensions
  const int ns[3] = {100, 100, 100}; //x, y, z gridpoint numbers
  double ds[3];  // distances between gridpoints
  double cycleErr;
  const int BCs[6] = {PoisFFT::PERIODIC, PoisFFT::PERIODIC,  //boundary conditions, all must be same type,
                      PoisFFT::PERIODIC, PoisFFT::PERIODIC,  //or else PoisFFT has undefined behaviour
                      PoisFFT::PERIODIC, PoisFFT::PERIODIC};
  int i;
  for (i = 0; i<3; i++){ // set the grid, depends on the boundary conditions
    ds[i] = Ls[i] / ns[i];
  }
  double *arr = new double[ns[0]*ns[1]*ns[2]]; // allocate the arrays contiguously, you can use any other class
  double *eps = new double[ns[0]*ns[1]*ns[2]]; //RELATIVE permittivity.
  double *chi = new double[ns[0]*ns[1]*ns[2]]; //electron affinity or WF
  double *RHS = new double[ns[0]*ns[1]*ns[2]]; // from which you can get a pointer to contiguous buffer, will contain rho.
  double *correction = new double[ns[0]*ns[1]*ns[2]]; // correction used to update RHS
  std::pair<int,double> *electrodemap = new std::pair<int,double>[ns[0]*ns[1]*ns[2]]; //stores electrode surface info and potentials.
  int cycleCount = 0;
  ProfilerStart("./profileresult.out"); //using google performance tools
  const std::clock_t begin_time = std::clock();
  init_eps(ns, ds, Ls, eps); // set permittivity
  init_rhs(ns, ds, Ls, chi, eps, RHS); // set rho and apply relative permittivity, also save electron affinity for bulk.
  // Electrodes elec1(0.02, 0.04, 0.02, 0.04, 0.03, 0.06, 0, WF_GOLD);
  // Electrodes elec2(0.02, 0.04, 0.06, 0.08, 0.03, 0.06, 0, WF_GOLD);
  // Electrodes elec3(0.06, 0.08, 0.02, 0.04, 0.03, 0.06, 10, WF_GOLD);
  // Electrodes elec4(0.06, 0.08, 0.06, 0.08, 0.03, 0.06, 20, WF_GOLD);
  Electrodes elec1(0.2, 0.4, 0.2, 0.4, 0.3, 0.6, 0, WF_GOLD);
  Electrodes elec2(0.2, 0.4, 0.6, 0.8, 0.3, 0.6, 0, WF_GOLD);
  Electrodes elec3(0.6, 0.8, 0.2, 0.4, 0.3, 0.6, 10, WF_GOLD);
  Electrodes elec4(0.6, 0.8, 0.6, 0.8, 0.3, 0.6, 10, WF_GOLD);
  // Electrodes elec1(2.0, 4.0, 2.0, 4.0, 3.0, 6.0, 5, WF_GOLD);
  // Electrodes elec2(2.0, 4.0, 6.0, 8.0, 3.0, 6.0, 10, WF_GOLD);
  // Electrodes elec3(6.0, 8.0, 2.0, 4.0, 3.0, 6.0, 20, WF_GOLD);
  // Electrodes elec4(6.0, 8.0, 6.0, 8.0, 3.0, 6.0, 30, WF_GOLD);
  // Electrodes elec1(0.0, 1.0, 0.0, 0.1, 0.3, 0.6, 0, WF_GOLD);
  // Electrodes elec2(0.0, 1.0, 0.9, 1.0, 0.3, 0.6, 0, WF_GOLD);
  // Electrodes elec3(0.2, 0.4, 0.4, 0.6, 0.3, 0.6, 10, WF_GOLD);
  // Electrodes elec4(0.6, 0.8, 0.4, 0.6, 0.3, 0.6, 10, WF_GOLD);
  // Electrodes elec1(0.6, 0.7, 0.1, 0.5, 0.3, 0.6, 10, WF_GOLD);
  elec1.draw(ns, ds, Ls, RHS, electrodemap, chi); //separately call draw for each electrode.
  elec2.draw(ns, ds, Ls, RHS, electrodemap, chi);
  elec3.draw(ns, ds, Ls, RHS, electrodemap, chi);
  elec4.draw(ns, ds, Ls, RHS, electrodemap, chi);
  PoisFFT::Solver<3, double> S(ns, Ls, BCs); //   create solver object, 3 dimensions, double precision
  std::cout << "Beginning solver" << std::endl;
  do{
    S.execute(arr, RHS); //run the solver, can be run many times for different right-hand side
    cycleErr = check_error(ns, ds, arr, correction, electrodemap, RHS, BCs);
    apply_correction(ns, ds, RHS, correction, electrodemap);
    std::cout << "On cycle " << cycleCount << " with error " << cycleErr << "." << std::endl;
    cycleCount++;
  }while(cycleErr > MAXERROR);
  std::cout << "Finished on cycle " << cycleCount << " with error " << cycleErr << std::endl;
  save_file2D(ns, ds, Ls, RHS, RHOFILE);
  save_file2D(ns, ds, Ls, arr, CORRECTIONFILE);
  save_file2D(ns, ds, Ls, arr, FILENAME); //solution is in arr
  ProfilerStop();
  std::cout << "Time elapsed: " << float(clock()-begin_time)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "Ending, deleting variables" << std::endl;
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
      RHS[i] += correction[i];
    }
  }
}

double check_error(const int ns[3], const double ds[3], double *arr, double *correction, std::pair<int,double> *electrodemap, double *RHS, const int BCs[6]){
  double err = 0;
  double d2Vtarget;
  double d2Vapprox;
  double xplust, xminust, yplust, yminust, zplust, zminust;
  double xplusa, xminusa, yplusa, yminusa, zplusa, zminusa;
  for(int i = 0; i < ns[0]; i++){
    for(int j = 0; j < ns[1]; j++){
      for(int k = 0; k < ns[2]; k++){
        if(electrodemap[IND(i,j,k)].first == true){ //only check error at electrodes.
          //for each side that is an electrode, add the electrode potential to correction.
          //for each side that is bulk, add its potential adjusted by the workfunction and electron affinity to correction
          // xminust = electrodemap[IND(i-1,j,k)].second;
          // yplust = electrodemap[IND(i,j+1,k)].second;
          // yminust = electrodemap[IND(i,j-1,k)].second;
          // zplust = electrodemap[IND(i,j,k+1)].second;
          // zminust = electrodemap[IND(i,j,k-1)].second;
          // xplusa = arr[IND(i+1,j,k)];
          // xminusa = arr[IND(i-1,j,k)];
          // yplusa = arr[IND(i,j+1,k)];
          // yminusa = arr[IND(i,j-1,k)];
          // zplusa = arr[IND(i,j,k+1)];
          // zminusa = arr[IND(i,j,k-1)];
          // if( BCs[0] == PoisFFT::PERIODIC ){ //all boundary conditions must be the same, all are PERIODIC. Wrap-around.
          //   if(i == ns[0]-1){
          //     xplust = electrodemap[IND(0,j,k)].second;
          //     xplusa = arr[IND(0,j,k)];
          //     // std::cout << "xmax" << std::endl;
          //   } else if(i == 0){
          //     xminust = electrodemap[IND(ns[0]-1,j,k)].second;
          //     xminusa = arr[IND(ns[0]-1,j,k)];
          //     // std::cout << "xmin" << std::endl;
          //   }
          //   if(j == ns[1]-1){
          //     yplust = electrodemap[IND(i,0,k)].second;
          //     yplusa = arr[IND(i,0,k)];
          //     // std::cout << "ymax" << std::endl;
          //   } else if(j == 0){
          //     yminust = electrodemap[IND(i,ns[1]-1,k)].second;
          //     yminusa = arr[IND(i,ns[1]-1,k)];
          //     // std::cout << "ymin" << std::endl;
          //   }
          //   if(k == ns[2]-1){
          //     zplust = electrodemap[IND(i,j,0)].second;
          //     zplusa = arr[IND(i,j,0)];
          //     // std::cout << "zmax" << std::endl;
          //   } else if(k == 0){
          //     zminust = electrodemap[IND(i,j,ns[2]-1)].second;
          //     zminusa = arr[IND(i,j,ns[2]-1)];
          //     // std::cout << "zmin" << std::endl;
          //   }
          // }
          d2Vtarget = (-6.0*electrodemap[IND(i,j,k)].second+electrodemap[IND(i-1,j,k)].second+electrodemap[IND(i+1,j,k)].second+
                      electrodemap[IND(i,j-1,k)].second+electrodemap[IND(i,j+1,k)].second+
                      electrodemap[IND(i,j,k-1)].second+electrodemap[IND(i,j,k+1)].second)/ds[0]/ds[0];
          // d2Vtarget = (-6.0*electrodemap[IND(i,j,k)].second+xplust+xminust+yplust+yminust+zplust+zminust)/ds[0]/ds[0];
          d2Vapprox = (-6.0*arr[IND(i,j,k)]+arr[IND(i-1,j,k)]+arr[IND(i+1,j,k)]+arr[IND(i,j-1,k)]+
                      arr[IND(i,j+1,k)]+arr[IND(i,j,k-1)]+arr[IND(i,j,k+1)])/ds[0]/ds[0];
          // d2Vapprox = (-6.0*arr[IND(i,j,k)]+xplusa+xminusa+yplusa+yminusa+zplusa+zminusa)/ds[0]/ds[0];
          correction[IND(i,j,k)] = d2Vtarget - d2Vapprox;
          err = std::max(err, fabs(correction[IND(i,j,k)]/d2Vtarget)); //get largest error value.
        }
      }
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
        " " << j*ds[1] << " " << arr[IND(i,j,k)] << std::endl;
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
          a[IND(i,j,k)] = 1;  //Free space
        } else{
          a[IND(i,j,k)] = 1e3; //Si permittivity
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
        a[IND(i,j,k)] = -1e10*Q_E/EPS0/eps[IND(i,j,k)]; //Set default set to bulk volume charge density.
        chi[IND(i,j,k)] = CHI_SI;
      }
    }
  }
  std::cout << "Finished RHS initialisation" << std::endl;
}
