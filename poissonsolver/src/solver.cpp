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
#include <fstream>
#include <cmath>
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

void Solver::set_val( double &val, double a ){
  val = a;
}

void Solver::init_val( std::vector<double> &vec, double val ){
  vec.resize(N[0]*N[1]*N[2]);
  std::fill(vec.begin(), vec.end(), val);
}

void Solver::init_rho( void ){
  rho.resize(N[0]*N[1]*N[2]);
  double x, y, z;
  for( int i = 0; i < N[0]; i++){
    x = i*L[0]/N[0];
    for( int j = 0; j < N[1]; j++){
      y = j*L[1]/N[1];
      for( int k = 0; k < N[2]; k++){
        z = k*L[2]/N[2];
        rho[i*N[1]*N[2] + j*N[2] + k] = Q_E*sin(x*2*PI/L[0]); //set rho with plane wave in x direction (one full wave)
/*
        if( (i+j+k)%2 == 0 ){
          //set periodic lattice (each occupied particle has unoccupied adjacent sites)
          rho[i*N[1]*N[2] + j*N[2] + k] = Q_E;
        }else{
          rho[i*N[1]*N[2] + j*N[2] + k] = 0;
        }
*/
      }
    }
  }
}

void Solver::init_eps( void ){
  eps.resize(N[0]*N[1]*N[2]);
  double x, y, z;
  unsigned int mask = 0;
  for( int i = 0; i < N[0]; i++){
    x = i*L[0]/N[0];
    for( int j = 0; j < N[1]; j++){
      y = j*L[1]/N[1];
      for( int k = 0; k < N[2]; k++){
        z = k*L[2]/N[2];
//        eps[i*N[1]*N[2] + j*N[2] + k] = 1;
        if ( x <= L[0]/2 ){
          eps[i*N[1]*N[2] + j*N[2] + k] = 1;
        } else {
          eps[i*N[1]*N[2] + j*N[2] + k] = 100;
        }
      }
    }
  }
}

void Solver::solve( void ){
  if( solvemethod == SOR ){
    poisson3DSOR();
  } else if( solvemethod == JACOBI ){
    poisson3DJacobi();
  } else if( solvemethod == GAUSS_SEIDEL ){
    poisson3DGaussSeidel();
  } else if( solvemethod == SOR_GEN ){
    poisson3DSOR_gen();
  }
}

void Solver::set_BCs (double Vx0, double VxL, double Vy0, double VyL, double Vz0, double VzL){
  int i = 0;
  int j = 0;
  int k = 0;
  for (j = 0; j<N[1]; j++){  //x = 0
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = Vx0;
      V[i*N[1]*N[2] + j*N[2] + k] = Vx0;
    }
  }
  i = N[0] - 1;  //x = Lx
  for (j = 0; j<N[1]; j++){
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = VxL;
      V[i*N[1]*N[2] + j*N[2] + k] = VxL;
    }
  }
  j = 0;  //y = 0
  for (i = 0; i<N[0]; i++){
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = Vy0;
      V[i*N[1]*N[2] + j*N[2] + k] = Vy0;
    }
  }
  j = N[1]-1;  //y = Ly
  for (i = 0; i<N[0]; i++){
    for (k = 0; k<N[2]; k++){
      rho[i*N[1]*N[2] + j*N[2] + k] = VyL;
      V[i*N[1]*N[2] + j*N[2] + k] = VyL;
    }
  }
  k = 0;  //z = 0
  for (i = 0; i<N[0]; i++){
    for (j = 0; j<N[1]; j++){
      rho[i*N[1]*N[2] + j*N[2] + k] = Vz0;
      V[i*N[1]*N[2] + j*N[2] + k] = Vz0;
    }
  }
  k = 0;  //z = Lz
  for (i = 0; i<N[0]; i++){
    for (j = 0; j<N[1]; j++){
      rho[i*N[1]*N[2] + j*N[2] + k] = VzL;
      V[i*N[1]*N[2] + j*N[2] + k] = VzL;
    }
  }
}

void Solver::write( std::vector<double> &vals, std::string filename ){
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping data to " << filename << std::endl;
  const std::vector<double> incSpacing = {L[0]/N[0], L[1]/N[1], L[2]/N[2]};
  for (int i = 0; i < N[0]; i++){
    for (int j = 0; j < N[1]; j++){
      for (int k = 0; k < N[2]; k++){
        //save data as x y z V
        outfile << std::setprecision(5) << std::scientific << i * incSpacing[0] << " " << j * incSpacing[1] <<
                " " << k * incSpacing[2] << " " << vals[i*N[1]*N[2] + j*N[2] + k] << std::endl;
      }
    }
    outfile << std::endl;
  }
  std::cout << "Ending." << std::endl;
}

void Solver::write_2D( std::vector<double> &vals, std::string filename ){
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping data to " << filename << std::endl;
  const std::vector<double> incSpacing = {L[0]/N[0], L[1]/N[1], L[2]/N[2]};
//  const int k = N[2]/2;
  const int k = N[2]/2;
  for (int i = 0; i < N[0]; i++){
    for (int j = 0; j < N[1]; j++){
        outfile << std::setprecision(5) << std::scientific << i * incSpacing[0] << " " << j * incSpacing[1] <<
                " " << vals[i*N[1]*N[2] + j*N[2] + k] << std::endl;
    }
    outfile << std::endl;
  }
  std::cout << "Ending." << std::endl;
}

std::vector<double> Solver::get_a( std::vector<double> &eps, const int &ind){
  std::vector<double> a(7);
  a[0] = (eps[ind] + eps[ind - 1] + eps[ind - N[2]] + eps[ind - N[2] - 1] +
         eps[ind - N[1]*N[2]] + eps[ind - N[1]*N[2] - 1] +
         eps[ind - N[1]*N[2] - N[2]] + eps[ind - N[1]*N[2] - N[2] - 1])/8;
  //i
  a[1] = (eps[ind] + eps[ind - 1] + eps[ind - N[2]] + eps[ind - N[2] - 1])/4;
  //j
  a[2] = (eps[ind] + eps[ind - 1] + eps[ind - N[1]*N[2]] + eps[ind - N[1]*N[2] - 1])/4;
    //k
  a[3] = (eps[ind] + eps[ind - N[2]] + eps[ind - N[1]*N[2]] + eps[ind - N[1]*N[2] - N[2]])/4;
    //i - 1
  a[4] = (eps[ind - N[1]*N[2]] + eps[ind - N[1]*N[2] - 1] +
         eps[ind - N[1]*N[2] - N[2]] + eps[ind - N[1]*N[2] - N[2] - 1])/4;
    //j - 1
  a[5] = (eps[ind - N[2]] + eps[ind - N[2] - 1] + eps[ind - N[1]*N[2] - N[2]] + eps[ind - N[1]*N[2] - N[2] - 1])/4;
  //k - 1
  a[6] = (eps[ind - 1] + eps[ind - N[2] - 1] + eps[ind - N[1]*N[2] - 1] + eps[ind - N[1]*N[2] - N[2] - 1])/4;
  return a;
}

void Solver::check_eps ( std::vector<double> &eps, std::vector<bool> &isChangingeps){
  int ind = N[1]*N[2] + N[2] + 1;
  for ( int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
    for ( int j = 1; j < N[1]-1; j++){
      for (int k = 1; k < N[2]-1; k++){
        if( (eps[ind] != eps[ind+1]) || (eps[ind] != eps[ind-1]) || (eps[ind] != eps[ind+N[2]]) ||
          (eps[ind] != eps[ind-N[2]]) || (eps[ind] != eps[ind+N[1]*N[2]]) || (eps[ind] != eps[ind-N[1]*N[2]]) ){
          //eps changing on a boundary
          isChangingeps[ind] = true;
        } else {
          //eps locally static
          isChangingeps[ind] = false;
        }
        ++ind;
      }
    ind += 2;
    }
  ind += 2*N[2];
  }
}

//based on http://www.eng.utah.edu/~cfurse/ece6340/LECTURE/FDFD/Numerical%20Poisson.pdf
void Solver::poisson3DSOR_gen( void ){
  double Vold; //needed to calculate error between new and old values
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
  const std::vector<double> overrelax = {1.85/6, -0.85}; //overrelaxation parameter for SOR
  std::vector<double> a(7);
  std::vector<bool> isChangingeps(N[0]*N[1]*N[2]);
  check_eps( eps, isChangingeps);
  std::cout << "DOING SOR_GEN" << std::endl;
  //obtain additional spacing value
  std::cout << "Iterating..." << std::endl;
  const std::clock_t begin_time = std::clock();
  unsigned int ind;
  do{
      ind = N[1]*N[2] + N[2] + 1;
      currError = 0;  //reset error for every run
      for ( int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
        for ( int j = 1; j < N[1]-1; j++){
          for (int k = 1; k < N[2]-1; k++){
            Vold = V[ind]; //Save for error comparison, can do in outer do{} loop to speed up.
            if( isChangingeps[ind] == true ){
              //there is a difference in permittivity, get a values
              a = get_a(eps, ind);
              V[ind] = overrelax[1]*V[ind] + overrelax[0]*(
              a[4]*V[ind-N[1]*N[2]] + a[1]*V[ind+N[1]*N[2]] +
              a[5]*V[ind-N[2]] + a[2]*V[ind+N[2]] +
              a[6]*V[ind-1] + a[3]*V[ind+1] + rho[ind]*h2/EPS0)/a[0]; //calculate new potential
            } else { //no difference in permittivity, do not calculate new a values
              V[ind] = overrelax[1]*V[ind] + overrelax[0]*(
              V[ind-1*N[1]*N[2]] + V[ind+1*N[1]*N[2]] +
              V[ind-1*N[2]] + V[ind+1*N[2]] +
              V[ind-1] + V[ind+1] + rho[ind]*h2/EPS0/eps[ind]); //calculate new potential
            }
            if ( fabs((V[ind] - Vold)/ V[ind]) > currError){ //capture worst case error
                currError = fabs((V[ind] - Vold)/V[ind]);
            }
            ++ind;
          }
        ind += 2;
        }
      ind += 2*N[2];
      }
      cycleCount++;
      if (cycleCount%50 == 0){
        std::cout << "On iteration " << cycleCount << " with " << currError << "% error." << std::endl;
      }
  }while( currError > MAXERROR );
  std::cout << "Finished in " << cycleCount << " iterations." << std::endl;
  std::cout << "Time elapsed: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << std::endl;
}

void Solver::poisson3DSOR( void ){
  double Vold;
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
  const std::vector<double> overrelax = {1.85/6, -0.85}; //overrelaxation parameter for SOR
  std::cout << "DOING SOR" << std::endl;
  //obtain additional spacing value
  std::cout << "Iterating..." << std::endl;
  const std::clock_t begin_time = std::clock();
  unsigned int ind;
  do{
      ind = N[1]*N[2] + N[2] + 1;
      //Vold = V; //Save for error comparison, save all before entering loop, so loop can focus on computation
      currError = 0;  //reset error for every run
      for ( int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
        for ( int j = 1; j < N[1]-1; j++){
          for (int k = 1; k < N[2]-1; k++){
            Vold = V[ind];
            V[ind] = overrelax[1]*V[ind] + overrelax[0]*(
            V[ind-1*N[1]*N[2]] + V[ind+1*N[1]*N[2]] +
            V[ind-1*N[2]] + V[ind+1*N[2]] +
            V[ind-1] + V[ind+1] + rho[ind]*h2/EPS0); //calculate new potential
            if ( fabs((V[ind] - Vold)/ V[ind]) > currError){ //capture worst case error
                currError = fabs((V[ind] - Vold)/V[ind]);
            }
          ++ind;
          }
        ind+= 2;
        }
      ind += 2*N[2];
      }
      cycleCount++;
      if (cycleCount%50 == 0){
        std::cout << "On iteration " << cycleCount << " with " << currError << "% error." << std::endl;
      }
  }while( currError > MAXERROR );
  std::cout << "Finished in " << cycleCount << " iterations." << std::endl;
  std::cout << "Time elapsed: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << std::endl;
}

void Solver::poisson3DJacobi( void ){
  std::vector<double> Vold; //needed to calculate error between new and old values
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
  Vold.resize(V.size());
  std::cout << "DOING JACOBI" << std::endl;
  //obtain additional spacing value
  std::cout << "Iterating..." << std::endl;
  const std::clock_t begin_time = std::clock();
  unsigned int ind;
  do{
      currError = 0;  //reset error for every run
      Vold = V; //Save for error comparison
      for ( int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
        for ( int j = 1; j < N[1]-1; j++){
          for (int k = 1; k < N[2]-1; k++){
            ind = i*N[1]*N[2] + j*N[2] + k;
            V[ind] = (Vold[ind-1*N[1]*N[2]] + Vold[ind+1*N[1]*N[2]] +
                     Vold[ind-1*N[2]] + Vold[ind+1*N[2]] +
                     Vold[ind-1] + Vold[ind+1] + rho[ind]*h2/EPS0)/6; //calculate new potential
            if ( fabs((V[ind] - Vold[ind])/ V[ind]) > currError){ //capture worst case error
                currError = fabs((V[ind] - Vold[ind])/V[ind]);
            }
          }
        }
      }
      cycleCount++;
      if (cycleCount%50 == 0){
        std::cout << "On iteration " << cycleCount << " with " << currError << "% error." << std::endl;
      }
  }while( currError > MAXERROR );
  std::cout << "Finished in " << cycleCount << " iterations." << std::endl;
  std::cout << "Time elapsed: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << std::endl;
}

void Solver::poisson3DGaussSeidel ( void ){
  double Vold; //needed to calculate error between new and old values
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
  std::cout << "DOING GAUSS SEIDEL" << std::endl;
  //obtain additional spacing value
  std::cout << "Iterating..." << std::endl;
  const std::clock_t begin_time = std::clock();
  unsigned int ind;
  do{
      currError = 0;  //reset error for every run
      for ( int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
        for ( int j = 1; j < N[1]-1; j++){
          for (int k = 1; k < N[2]-1; k++){
            ind = i*N[1]*N[2] + j*N[2] + k;
            Vold = V[ind]; //Save for error comparison
            V[ind] = (V[ind-1*N[1]*N[2]] + V[ind+1*N[1]*N[2]] +
                     V[ind-1*N[2]] + V[ind+1*N[2]] +
                     V[ind-1] + V[ind+1] + rho[ind]*h2/EPS0)/6; //calculate new potential
            if ( fabs((V[ind] - Vold)/ V[ind]) > currError){ //capture worst case error
                currError = fabs((V[ind] - Vold)/V[ind]);
            }
          }
        }
      }
      cycleCount++;
      if (cycleCount%50 == 0){
        std::cout << "On iteration " << cycleCount << " with " << currError << "% error." << std::endl;
      }
  }while( currError > MAXERROR );
  std::cout << "Finished in " << cycleCount << " iterations." << std::endl;
  std::cout << "Time elapsed: " << float( clock () - begin_time ) /  CLOCKS_PER_SEC << " seconds" << std::endl;
}
