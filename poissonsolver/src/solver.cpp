//
// Created by nathan on 03/05/17.
//
//Compute solution to Poisson's equation assuming Dirichlet boundary condition.
//

#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <functional>
#include "solver.h"
#include <ctime>

#define MAXERROR 1e-5

//Solver class functions
void Solver::set_val( std::vector<int> &vec, int a, int b, int c ){
  vec.resize(3);
  vec[0] = a;
  vec[1] = b;
  vec[2] = c;
}

void Solver::set_val( std::vector<double> &vec, double a, double b, double c ){
  vec.resize(3);
  vec[0] = a;
  vec[1] = b;
  vec[2] = c;
}

void Solver::set_val( int &val, int a ){
  val = a;
}

void Solver::init_val( std::vector<double> &vec, double val ){
  vec.resize(N[0]*N[1]*N[2]);
  std::fill(vec.begin(), vec.end(), val);
}

void Solver::solve( void ){
  if( solvemethod == SOR ){
    poisson3DSOR();
  }
}

void Solver::set_BCs (double Vx0, double VxL, double Vy0, double VyL, double Vz0, double VzL){
  int i = 0;
  int j = 0;
  int k = 0;
  //x = 0
  for (j = 0; j<N[1]; j++){
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = Vx0;
      V[i*N[1]*N[2] + j*N[2] + k] = Vx0;
    }
  }
  //x = Lx
  i = N[0] - 1;
  for (j = 0; j<N[1]; j++){
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = VxL;
      V[i*N[1]*N[2] + j*N[2] + k] = VxL;
    }
  }
  //y = 0
  j = 0;
  for (i = 0; i<N[0]; i++){
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = Vy0;
      V[i*N[1]*N[2] + j*N[2] + k] = Vy0;
    }
  }
  //y = Ly
  j = N[1]-1;
  for (i = 0; i<N[0]; i++){
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = VyL;
      V[i*N[1]*N[2] + j*N[2] + k] = VyL;
    }
  }
  //z = 0
  k = 0;
  for (i = 0; i<N[0]; i++){
    for (j = 0; j<N[1]; j++){
      rho[i*N[1]*N[2] + j*N[2] + k] = Vz0;
      V[i*N[1]*N[2] + j*N[2] + k] = Vz0;
    }
  }
  //z = Lz
  k = 0;
  for (i = 0; i<N[0]; i++){
    for (j = 0; j<N[1]; j++){
      rho[i*N[1]*N[2] + j*N[2] + k] = VzL;
      V[i*N[1]*N[2] + j*N[2] + k] = VzL;
    }
  }
}
void Solver::poisson3DSOR( void ){
  /*
    //Parameters
      std::cout << "Setup for SOR" << std::endl;
  */

  //N - lattice points in x, y, z
  //L - sample lengths in x, y, z
  //V - vector for potential, sized correctly if init_val was called
  //rho - charge density, sized correctly if init_val was called
  //
  //MAXERROR defined at top
  //EPS0 defined at top
  //
  double Vold; //needed to calculate error between new and old values
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
  double overrelax = 1.85; //overrelaxation parameter for SOR
  std::cout << "DOING SOR" << std::endl;

  //obtain additional spacing value
  h2.push_back(1/(h2[0]+h2[1]+h2[2])/2);
  std::cout << "Iterating..." << std::endl;
  const std::clock_t begin_time = std::clock();
  do{
      currError = 0;  //reset error for every run
      for ( int i = 1; i < N[0]-1; i=i+1){ //for all x points except endpoints
        for ( int j = 1; j < N[1]-1; j=j+1){
          for (int k = 1; k < N[2]-1; k=k+1){
            Vold = V[i*N[1]*N[2] + j*N[2] + k]; //Save for error comparison
            V[i*N[1]*N[2] + j*N[2] + k] = (1-overrelax)*V[i*N[1]*N[2] + j*N[2] + k] + overrelax*(
              (V[(i-1)*N[1]*N[2] + j*N[2] + k] + V[(i+1)*N[1]*N[2] + j*N[2] + k])*h2[0] +
              (V[i*N[1]*N[2] + (j-1)*N[2] + k] + V[i*N[1]*N[2] + (j+1)*N[2] + k])*h2[1] +
              (V[i*N[1]*N[2] + j*N[2] + (k+1)] + V[i*N[1]*N[2] + j*N[2] + (k-1)])*h2[2] +
              rho[i*N[1]*N[2] + j*N[2] + k])*h2[3]; //calculate new potential
            if ( fabs((V[i*N[1]*N[2] + j*N[2] + k] - Vold)/V[i*N[1]*N[2] + j*N[2] + k]) > currError ){ //capture worst case error
                currError = fabs((V[i*N[1]*N[2] + j*N[2] + k] - Vold)/V[i*N[1]*N[2] + j*N[2] + k]);
            }
          }
        }
      }
      cycleCount++;
      if (cycleCount%50 == 0){
        std::cout << "On iteration " << cycleCount << " with " << currError << "% error." << std::endl;
      }
  }while( currError > MAXERROR );

  std::cout << "Finished" << std::endl;
  std::cout << "Time elapsed: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << std::endl;
}

//Reference poisson functions
void poisson1DJacobi(void) {
    //Parameters
    unsigned int N = 1000;           //Number of lattice points
    double sampleLength = 1;   //Physical length of sample in meters
    std::vector<double> V(N,0);     //Electric potential solution to Poisson's equation (Volts)
    std::vector<double> Vold(N,0);     //Vector to copy old values
    std::vector<double> rho(N,1.6e-19);   //Charge density, first and last are boundary conditions
    double maxError = 1e-5;       //set at 1% error allowed
    double currError;
    unsigned int cycleCount = 0;
    //Setup boundary condition
    std::cout << "Setup for Jacobi" << std::endl;
    rho[0] = 0;
    rho[N-1] = 0;
    //Prepare rho with EPS and spacing multiplication
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::divides<double>(), EPS0 ) );
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::multiplies<double>(),(sampleLength/(N-1)*(sampleLength/(N-1)))));
    //Apply coundary condition
    std::cout << "Applying boundary condition" << std::endl;
    V[0] = rho[0];
    V[N-1] = rho[N-1];
    //Computation loop
    std::cout << "Iterating..." << std::endl;
    do{
        currError = 0;  //reset error for every run
        Vold = V; //need to copy over all old values for Jacobi
        for( int i = 1; i < N-1; i++){ //for all lattice points except endpoints
            V[i] = (Vold[i-1] + Vold[i+1] + rho[i])/2.0; //calculate new potential
            if ( fabs((V[i] - Vold[i])/V[i]) > currError ){ //capture worst case error
                currError = fabs((V[i] - Vold[i])/V[i]);
            }
        }
        cycleCount++;
/*
        if( cycleCount%1000  == 0 ){
          std::cout << "On iteration number " << cycleCount << std::endl;
        }
*/
    }while( currError > maxError );

    std::cout << "Finished in " << cycleCount << " iterations within " <<
              maxError*100 << "% error." << std::endl;
    int i = 0;
    while ( i < N ){
//        std::cout << "V[" << i << "] = " << std::setprecision(20) << V[i] << "V" <<std::endl;
//        std::cout << i << "," << V[i] <<std::endl; //csv-friendly
        i++;
    }
}

void poisson1DGaussSeidel(void) {
    //Parameters
    unsigned int N = 1000;           //Number of lattice points
    double sampleLength = 1;       //Physical length of sample in meters
    std::vector<double> V(N,0);    //Electric potential solution to Poisson's equation (Volts)
    std::vector<double> rho(N,1.6e-19);   //Charge density, first and last are boundary conditions
    double Vold;                   //Vold only requires one element to calculate error
    double maxError = 1e-5;        //set at 1% error allowed
    double currError;
    unsigned int cycleCount = 0;
    //Setup
    std::cout << "Setup for Gauss-Seidel" << std::endl;
    rho[0] = 0;
    rho[N-1] = 0;
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::divides<double>(), EPS0 ) );
    std::transform( rho.begin()+1, rho.end()-1, //multiply by spacing squared
                    rho.begin()+1, std::bind2nd(std::multiplies<double>(),(sampleLength/(N-1)*(sampleLength/(N-1)))));
    V[0] = rho[0];
    V[N-1] = rho[N-1];
    //Computation loop
    std::cout << "Iterating..." << std::endl;
    do{
        currError = 0;  //reset error for every run
        for( int i = 1; i < N-1; i=i+1 ){ //for all lattice points except endpoints,
            Vold = V[i]; //Save for error comparison
            V[i] = (V[i-1] + V[i+1] + rho[i])/2.0; //calculate new potential
            if ( fabs((V[i] - Vold)/V[i]) > currError ){ //capture worst case error
                currError = fabs((V[i] - Vold)/V[i]);
            }
        }
        cycleCount++;
/*
        if( cycleCount%1000  == 0 ){
          std::cout << "On iteration number " << cycleCount << std::endl;
        }
*/
    }while( currError > maxError );

    std::cout << "Finished in " << cycleCount << " iterations within " <<
              maxError*100 << "% error." << std::endl;
    int i = 0;
    while ( i < N ){
//        std::cout << "V[" << i << "] = " << std::setprecision(20) << V[i] << "V" <<std::endl;
//        std::cout << i << "," << V[i] <<std::endl; //csv-friendly
        i++;
    }
}

void poisson1DSOR(void) {
    //Parameters
    unsigned int N = 1000;           //Number of lattice points
    double sampleLength = 1;       //Physical length of sample in meters
    std::vector<double> V(N,0);    //Electric potential solution to Poisson's equation (Volts)
    std::vector<double> rho(N,1.6e-19);   //Charge density, first and last are boundary conditions
    double Vold;                   //Vold only requires one element to calculate error
    double maxError = 1e-5;        //set at 1% error allowed
    double overrelax = 1.85;        //overrelaxation parameter for SOR method
    double currError;
    unsigned int cycleCount = 0;
    //Setup
    std::cout << "Setup for SOR" << std::endl;
    rho[0] = 0;
    rho[N-1] = 0;
    std::transform( rho.begin()+1, rho.end()-1,
                    rho.begin()+1, std::bind2nd(std::divides<double>(), EPS0 ) );
    std::transform( rho.begin()+1, rho.end()-1, //multiply by spacing squared
                    rho.begin()+1, std::bind2nd(std::multiplies<double>(),(sampleLength/(N-1)*(sampleLength/(N-1)))));
    V[0] = rho[0];
    V[N-1] = rho[N-1];
    //Computation loop
    std::cout << "Iterating..." << std::endl;
    do{
        currError = 0;  //reset error for every run
        for( int i = 1; i < N-1; i=i+1 ){ //for all lattice points except endpoints,
            Vold = V[i]; //Save for error comparison
            V[i] = (1-overrelax)*V[i] + overrelax*(V[i-1] + V[i+1] + rho[i])/2.0; //calculate new potential
            if ( fabs((V[i] - Vold)/V[i]) > currError ){ //capture worst case error
                currError = fabs((V[i] - Vold)/V[i]);
            }
        }
        cycleCount++;
/*
        if( cycleCount%1000  == 0 ){
          std::cout << "On iteration number " << cycleCount << std::endl;
        }
*/
    }while( currError > maxError );

    std::cout << "Finished in " << cycleCount << " iterations within " <<
              maxError*100 << "% error." << std::endl;
    int i = 0;
    while ( i < N ){
//        std::cout << "V[" << i << "] = " << std::setprecision(20) << V[i] << "V" <<std::endl;
//        std::cout << i << "," << V[i] <<std::endl; //csv-friendly
        i++;
    }
}
