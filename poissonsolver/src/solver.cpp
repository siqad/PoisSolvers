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

#define EPS0 8.85418782e-12
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
    poisson3DSOR( rho );
  }
}

void Solver::poisson3DSOR( std::vector<double> rho){

  /*
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
  */

  //N - lattice points in x, y, z
  //L - sample lengths in x, y, z
  //V - vector for potential, sized correctly if init_V was called
  //
  //MAXERROR defined at top
  //rho - charge density, provided as argument
  //EPS0 defined at top
  //
  double Vold; //needed to calculate error between new and old values
  double currError = 0;
  unsigned int cycleCount = 0;
  double overrelax = 1.85;

  std::cout << "DOING SOR" << std::endl;




/*
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

          if( cycleCount%1000  == 0 ){
            std::cout << "On iteration number " << cycleCount << std::endl;
          }

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
*/






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
