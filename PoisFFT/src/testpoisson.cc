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
*/

#include <cmath>
#include <iostream>
#include "poisfft.h"
#include <stdlib.h>
#include <iomanip>
#include <fstream>

const double pi = 3.14159265358979323846;

#define IND(i,j,k) (i)*(ns[1]*ns[2])+(j)*(ns[2])+k
#define FILENAME "outfile.txt"

void init_rhs(const int ns[3], const double ds[3], const double Ls[3], double* a){
  int i,j,k;
  std::cout << "Initialising RHS" << std::endl;
  for (i=0;i<ns[0];i++){
    for (j=0;j<ns[1];j++){
      for (k=0;k<ns[2];k++){
        double x = ds[0]*(i+0.5);
        double y = ds[1]*(j+0.5);
        double z = ds[2]*(k+0.5);
        //piece wise RHS
        //a[IND(i,j,k)] = sin(3*pi*x/Ls[0]) * sin(5*pi*y/Ls[1]) *sin(7*pi*z/Ls[2]);
        if (k < ns[2]/2){
        a[IND(i,j,k)] = x/Ls[0] * y/Ls[1] * z/Ls[2];;
        }
        else{
          a[IND(i,j,k)] = 0;
        }
      }
    }
  }
  std::cout << "Finished initialisation" << std::endl;
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
                   (- 1.0 / ( pow(3*pi/Ls[0], 2.0) +
                              pow(5*pi/Ls[1], 2.0) +
                              pow(7*pi/Ls[2], 2.0) ) *
                    sin(3*pi*x/Ls[0]) *
                    sin(5*pi*y/Ls[1]) *
                    sin(7*pi*z/Ls[2]))
                , 2.0);
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

int main(){
  std::cout << "Modified by Nathan Chiu. Code is offered as is, with no warranty." <<
    "See LICENCE.GPL for licence info." << std::endl;
  // testing file (csv)
  std::ofstream outfile;
  outfile.open("outfile.txt", std::ios_base::out | std::ios_base::trunc );
  // domain dimensions
  const double Ls[3] = {1.0, 2.0, 3.0}; //x, y, z

  // gridpoint numbers
  const int ns[3] = {20, 20, 20}; //x, y, z

  // distances between gridpoints
  double ds[3];

  //boundary conditions
/*
  const int BCs[6] = {PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG,
                      PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG,
                      PoisFFT::NEUMANN_STAG, PoisFFT::NEUMANN_STAG};
*/
  const int BCs[6] = {PoisFFT::DIRICHLET, PoisFFT::DIRICHLET,
                      PoisFFT::DIRICHLET, PoisFFT::DIRICHLET,
                      PoisFFT::DIRICHLET, PoisFFT::DIRICHLET};
  int i;
  // set the grid, depends on the boundary conditions
  for (i = 0; i<3; i++){
    ds[i] = Ls[i] / ns[i];
  }

  // allocate the arrays contiguously, you can use any other class
  // from which you can get a pointer to contiguous buffer
  double *arr = new double[ns[0]*ns[1]*ns[2]];
  double *RHS = new double[ns[0]*ns[1]*ns[2]];

  // set the right-hand side
  init_rhs(ns, ds, Ls, RHS);

  // create solver object, 3 dimensions, double precision
  PoisFFT::Solver<3, double> S(ns, Ls, BCs);

  //run the solver, can be run many times for different right-hand sides
  S.execute(arr, RHS);

  // check correctness (compares with known, exact solution)
  //check_solution(ns, ds, Ls, arr);

  //Dump to file. Maybe consider adding as a class method?
  std::cout << "Dumping to " << FILENAME << std::endl;
  for (int i = 0; i < ns[0]; i++){
    for (int j = 0; j < ns[1]; j++){
      for (int k = 0; k < ns[2]; k++) {
        //std::cout << arr[ j*ns[1]*ns[2] + k*ns[2] + l ] << "," << std::endl;
        //save data as x y z V
        outfile << std::setprecision(5) << std::fixed << i * ds[0] << " " << j * ds[1] <<
                " " << k * ds[2] << " " << arr[i * ns[1] * ns[2] + j * ns[2] + k] << std::endl;
      }
    }
    outfile << std::endl;
  }
  std::cout << "Ending, deleting variables" << std::endl;

  delete[] RHS;
  delete[] arr;
}
