//
// Created by nathan on 03/05/17.
// A BUNCH OF SOLVER CLASS METHODS
// Compute solution to Poisson's equation assuming Dirichlet boundary condition.
// Maybe seal off Jacobi and Gauss-Seidel and normal SOR methods in favor of SOR_GEN?

#include "solver.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <ctime>
#include <fstream>
#include <cmath>


//set_val, used to assign values to class attributes.
int* Solver::set_val( int a, int b, int c ){
  int * arr = new int[3];
  arr[0] = a;
  arr[1] = b;
  arr[2] = c;
  return arr;
}
double* Solver::set_val( double a, double b, double c ){
  double * arr = new double[3];
  arr[0] = a;
  arr[1] = b;
  arr[2] = c;
  return arr;
}
double* Solver::set_val( double a, double b){
  double * arr = new double[3];
  arr[0] = a;
  arr[1] = b;
  return arr;
}

//init_val, used for initializing arrays.
double * Solver::init_val( double val, double* prev ){ //fills the entire vector with val
  double * arr = new double[N[0]*N[1]*N[2]];
  for(int i = 0; i < N[0]*N[1]*N[2]; i++){
    arr[i] = val;
  }
  return arr;
}
std::pair<double,double>* Solver::init_val( double val, std::pair<double, double>* prev ){ //fills the entire vector with val
  std::pair<double,double> * arr = new std::pair<double, double>[N[0]*N[1]*N[2]];
  for(int i = 0; i < N[0]*N[1]*N[2]; i++){
    arr[i].first = val;
    arr[i].second = val;
  }
  return arr;
}

//del, used for deleting arrays crated with new
void Solver::del( double* arr ){
  delete[] arr;
}
void Solver::del( int* arr ){
  delete[] arr;
}
void Solver::del( std::pair<double, double>* arr ){
  delete[] arr;
}

//init_rho, initializes rho with values.
//TODO change to read from file or from rho object or pointer to rho object
void Solver::init_rho( void ){
  double x, y, z;
  std::cout << "Initialising rho..." << std::endl;
  for(int i = 0; i < N[0]; i++){
    x = i*L[0]/N[0];
    for(int j = 0; j < N[1]; j++){
      y = j*L[1]/N[1];
      for(int k = 0; k < N[2]; k++){
        z = k*L[2]/N[2]; //set rho with plane wave in x direction (one full wave)
          rho[i*N[1]*N[2]+j*N[2]+k] = 1e10*Q_E*sin(x*2*PI/L[0]);
      }
    }
  }
}

//init_eps, initializes eps with values.
//TODO change to read from file or from eps object or pointer to eps object
void Solver::init_eps( void ){
  std::cout << "Initialising eps..." << std::endl;
  eps = new double[N[0]*N[1]*N[2]];
  double x, y, z;
  for( int i = 0; i < N[0]; i++){
    x = i*L[0]/N[0];
    for( int j = 0; j < N[1]; j++){
      y = j*L[1]/N[1];
      for( int k = 0; k < N[2]; k++){
        z = k*L[2]/N[2];
        eps[i*N[1]*N[2]+j*N[2]+k] = 10*sin(x*2*PI/L[0])*sin(y*2*PI/L[1]) + 15;
      }
    }
  }
}

//calls the appropriate poisson solver.
void Solver::solve( void ){
  switch ( solvemethod ) {
    case SOR:           poisson3DSOR();
                        break;
    case JACOBI:        poisson3DJacobi();
                        break;
    case GAUSS_SEIDEL:  poisson3DGaussSeidel();
                        break;
    case SOR_GEN:       poisson3DSOR_gen();
  }
}

//sets boundary condition for both Dirichlet and Neumann types.
void Solver::set_BCs (double Vx0, double VxL, double Vy0, double VyL, double Vz0, double VzL){
  int i = 0;
  int j = 0;//RHO contains the boundary condition. setting V is for Dirichlet.
  int k = 0;//if Neumann BC, V will be overwritten in calc_Neumann()
  std::cout << "Applying boundary conditions..." << std::endl;
  for (j = 0; j < N[1]; j++){  //x = 0
    for (k = 0; k < N[2]; k++){
      rho[i*N[1]*N[2]+j*N[2]+k] = Vx0;
    }
  }
  i = N[0]-1;  //x = Lx
  for (j = 0; j < N[1]; j++){
    for (k = 0; k < N[2]; k++){
      rho[i*N[1]*N[2]+j*N[2]+k] = VxL;
      V[i*N[1]*N[2]+j*N[2]+k] = VxL;
    }
  }
  j = 0;  //y = 0
  for (i = 0; i < N[0]; i++){
    for (k = 0; k < N[2]; k++){
      rho[i*N[1]*N[2]+j*N[2]+k] = Vy0;
      V[i*N[1]*N[2]+j*N[2]+k] = Vy0;
    }
  }
  j = N[1]-1;  //y = Ly
  for (i = 0; i < N[0]; i++){
    for (k = 0; k < N[2]; k++){
      rho[i*N[1]*N[2]+j*N[2]+k] = VyL;
      V[i*N[1]*N[2]+j*N[2]+k] = VyL;
    }
  }
  k = 0;  //z = 0
  for (i = 0; i < N[0]; i++){
    for (j = 0; j < N[1]; j++){
      rho[i*N[1]*N[2]+j*N[2]+k] = Vz0;
      V[i*N[1]*N[2]+j*N[2]+k] = Vz0;
    }
  }
  k = 0;  //z = Lz
  for (i = 0; i < N[0]; i++){
    for (j = 0; j < N[1]; j++){
      rho[i*N[1]*N[2]+j*N[2]+k] = VzL;
      V[i*N[1]*N[2]+j*N[2]+k] = VzL;
    }
  }
}

//write 3D data to file
void Solver::write( std::vector<double> &vals, std::string filename ){
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping data to " << filename << std::endl;
  static const std::vector<double> incSpacing = {L[0]/N[0], L[1]/N[1], L[2]/N[2]};
  for (int i = 0; i < N[0]; i++){
    for (int j = 0; j < N[1]; j++){
      for (int k = 0; k < N[2]; k++){ //save data as x y z V
        outfile << std::setprecision(5) << std::scientific << i * incSpacing[0] << " " << j * incSpacing[1] <<
                " " << k * incSpacing[2] << " " << vals[i*N[1]*N[2] + j*N[2] + k] << std::endl;
      }
    }
    outfile << std::endl;
  }
  std::cout << "Ending." << std::endl;
}

//write 2D data to file
void Solver::write_2D( std::vector<double> &vals, std::string filename ){
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping data to " << filename << std::endl;
  static const std::vector<double> incSpacing = {L[0]/N[0], L[1]/N[1], L[2]/N[2]};
  static const int k = N[2]/2;
  for (int i = 0; i < N[0]; i++){
    for (int j = 0; j < N[1]; j++){
        outfile << std::setprecision(5) << std::scientific << i * incSpacing[0] << " " << j * incSpacing[1] <<
                " " << vals[i*N[1]*N[2] + j*N[2] + k] << std::endl;
    }
    outfile << std::endl;
  }
  std::cout << "Ending." << std::endl;
}
void Solver::write_2D( double* vals, std::string filename ){
  std::ofstream outfile;
  outfile.open(filename, std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping data to " << filename << std::endl;
  static const int k = N[2]/2;
  for (int i = 0; i < N[0]; i++){
    for (int j = 0; j < N[1]; j++){
        outfile << std::setprecision(5) << std::scientific << i*L[0]/N[0] << " " << j*L[1]/N[1] <<
                " " << vals[i*N[1]*N[2]+j*N[2]+k] << std::endl;
    }
    outfile << std::endl;
  }
  std::cout << "Ending." << std::endl;
}

//get new values for a, needed when permittivity changes locally.
void Solver::get_a( double* a, double* eps, int ind){
  a[0] = (eps[ind]+eps[ind-1]+eps[ind-N[2]]+eps[ind-N[2]-1]+eps[ind-N[1]*N[2]]+eps[ind-N[1]*N[2]-1]+
         eps[ind-N[1]*N[2]-N[2]]+eps[ind-N[1]*N[2]-N[2]-1])/8;
  a[1] = (eps[ind]+eps[ind-1]+eps[ind-N[2]]+eps[ind-N[2]-1])/4;//i
  a[2] = (eps[ind]+eps[ind-1]+eps[ind-N[1]*N[2]]+eps[ind-N[1]*N[2]-1])/4;//j
  a[3] = (eps[ind]+eps[ind-N[2]]+eps[ind-N[1]*N[2]]+eps[ind-N[1]*N[2]-N[2]])/4;//k
  a[4] = (eps[ind-N[1]*N[2]]+eps[ind-N[1]*N[2]-1]+eps[ind-N[1]*N[2]-N[2]]+eps[ind-N[1]*N[2]-N[2]-1])/4;//i - 1
  a[5] = (eps[ind-N[2]]+eps[ind-N[2]-1]+eps[ind-N[1]*N[2]-N[2]]+eps[ind-N[1]*N[2]-N[2]-1])/4;//j - 1
  a[6] = (eps[ind-1]+eps[ind-N[2]-1]+eps[ind-N[1]*N[2]-1]+eps[ind-N[1]*N[2]-N[2]-1])/4;//k - 1
}

//check whether permittivity changes locally, and remember for later
void Solver::check_eps (double *eps, bool* isChangingeps){
  int ind = N[1]*N[2]+N[2]+1;
  std::cout << "Checking eps..." << std::endl;
  for (int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
    for (int j = 1; j < N[1]-1; j++){
      for (int k = 1; k< N [2]-1; k++){
        if( (eps[ind] != eps[ind+1]) || (eps[ind] != eps[ind-1]) || (eps[ind] != eps[ind+N[2]]) ||
          (eps[ind] != eps[ind-N[2]]) || (eps[ind] != eps[ind+N[1]*N[2]]) || (eps[ind] != eps[ind-N[1]*N[2]]) ){
          isChangingeps[ind] = true; //eps changing on a boundary
        } else {
          isChangingeps[ind] = false; //eps changing on a boundary
        }
        ++ind;
      }
    ind += 2;
    }
  ind += 2*N[2];
  }
}

//check whether branch should do boundary or interior calculation
void Solver::check_exterior( bool* isExterior ){
  std::cout << "Checking exterior..." << std::endl;
  for (int i = 0; i < N[0]; i++){
    for (int j = 0; j < N[1]; j++){
      for (int k = 0; k < N[2]; k++){
        int ind = i*N[1]*N[2] + j*N[2] + k;
        if((i == 0) || (i == N[0]-1) || (j== N[1]-1) || (j == 0) || (k == 0) || (k == N[2]-1)){
          isExterior[ind] = true;
        } else {
          isExterior[ind] = false;
        }
      }
    }
  }
}

//check whether branch should do normal or ohmic calculation
void Solver::check_elec( bool* isBesideElec ){
  std::cout << "Checking electrodes..." << std::endl;
  for (int i = 0; i < N[0]; i++){
    for (int j = 0; j < N[1]; j++){
      for (int k = 0; k < N[2]; k++){
        int ind = i*N[1]*N[2] + j*N[2] + k;
        if( electrodemap[ind+1].first!=0 || electrodemap[ind-1].first!=0 ||
            electrodemap[ind+N[2]].first!=0 || electrodemap[ind-N[2]].first!=0 ||
            electrodemap[ind+N[1]*N[2]].first!=0 || electrodemap[ind-N[1]*N[2]].first!=0 ){
          isBesideElec[ind] = true;
        } else {
          isBesideElec[ind] = false;
        }
      }
    }
  }
}

//get value at boundary based on Neumann boundary condition
void Solver::calc_Neumann( int i, int j, int k){
  //For Neumann Boundary Conditions
    if( i == 0 ){ //x == 0;
      V[j*N[2] + k] = V[N[1]*N[2] + j*N[2] + k] - rho[j*N[2] + k]*h;
    } else if ( i == N[0]-1 ){ //x == Lx;
      V[(N[0]-1)*N[1]*N[2] + j*N[2] + k] = V[((N[0]-1)-1)*N[1]*N[2] + j*N[2] + k] - rho[(N[0]-1)*N[1]*N[2] + j*N[2] + k]*h;
    } else if ( j == 0 ) { //y == 0;
      V[i*N[1]*N[2] + k] = V[i*N[1]*N[2] + N[2] + k] - rho[i*N[1]*N[2] + k]*h;
    } else if ( j == N[1]-1 ){ //y == Ly;
      V[i*N[1]*N[2] + (N[1]-1)*N[2] + k] = V[i*N[1]*N[2] + ((N[1]-1)-1)*N[2] + k] - rho[i*N[1]*N[2] + (N[1]-1)*N[2] + k]*h;
    } else if ( k == 0) { //z == 0;
      V[i*N[1]*N[2] + j*N[2]] = V[i*N[1]*N[2] + j*N[2] + 1] - rho[i*N[1]*N[2] + j*N[2]]*h;
    } else { //z == Lz;
      V[i*N[1]*N[2] + j*N[2] + (N[2]-1)] = V[i*N[1]*N[2] + j*N[2] + ((N[2]-1)-1)] - rho[i*N[1]*N[2] + j*N[2] + (N[2]-1)]*h;
    }
}

void Solver::create_a( double** a){
  unsigned long int ind;
  std::cout << "Initialising a[0.." << N[0]*N[1]*N[2] << "][0..7]" << std::endl;
  for (int i = 0; i < N[0]*N[1]*N[2]; ++i){
    a[i] = new double[7]; //7 a values per point.
  }
  for ( int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
    for ( int j = 1; j < N[1]-1; j++){
      for (int k = 1; k < N[2]-1; k++){
        ind = i*N[1]*N[2] + j*N[2] + k;
        get_a(a[ind], eps, ind);
      }
    }
  }
}

//uses Generalised Successive Over Relaxation method to solve poisson's equation.
//based on http://www.eng.utah.edu/~cfurse/ece6340/LECTURE/FDFD/Numerical%20Poisson.pdf
void Solver::poisson3DSOR_gen( void ){
  double Vold; //needed to calculate error between new and old values
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
  double * overrelax = set_val(1.85/6, -0.85);
//  double * a = new double[7];
  double **a = new double*[N[0]*N[1]*N[2]]; //array of pointers to doubles.
  bool* isBesideElec = new bool[N[0]*N[1]*N[2]];
  bool* isExterior = new bool[N[0]*N[1]*N[2]]; //precompute whether or not exterior point.
  bool* isChangingeps = new bool[N[0]*N[1]*N[2]];
  check_eps( eps, isChangingeps );
  check_exterior( isExterior );
  check_elec( isBesideElec );
  create_a( a );
  std::cout << "DOING SOR_GEN" << std::endl;
  std::cout << "Iterating..." << std::endl;
  const std::clock_t begin_time = std::clock();
  unsigned long int ind;
  do{
      currError = 0;  //reset error for every run
      for ( int i = 0; i < N[0]; i++){ //for all x points except endpoints
        for ( int j = 0; j < N[1]; j++){
          for (int k = 0; k < N[2]; k++){
            ind = i*N[1]*N[2] + j*N[2] + k;
            //directly on face, do boundary first.
            if ((boundarytype == NEUMANN)&&(isExterior[ind]==true)){
                calc_Neumann(i, j, k); //Dirichlet boundary handled by set_BCs() in main
            } else if (isExterior[ind] == false){ //interior point, do normal.
              if( electrodemap[ind].first == 0){ //current cell is NOT electrode, perform calculation
                Vold = V[ind]; //Save for error comparison
                if(isBesideElec[ind] == true){
                //Current cell is NEXT TO an electrode, ignore other effects and use workfunction ohmic contact calculation
                //Rectangular electrodes, only one side of bulk can interface with electrode.
                //Work function and electrode voltage = 0 for non-electrode sites. Add all sides, since only one of them
                //is non-zero
                  V[ind] = (electrodemap[ind+1].second + electrodemap[ind-1].second +
                           electrodemap[ind+N[2]].second + electrodemap[ind-N[2]].second +
                           electrodemap[ind+N[1]*N[2]].second + electrodemap[ind-N[1]*N[2]].second) -
                           ((electrodemap[ind+1].first + electrodemap[ind-1].first +
                           electrodemap[ind+N[2]].first + electrodemap[ind-N[2]].first +
                           electrodemap[ind+N[1]*N[2]].first + electrodemap[ind-N[1]*N[2]].first) - CHI_SI);
                } else { //not directly beside an electrode, perform normal calculation.
                  if( isChangingeps[ind] == true ){ //check if at a permittivity boundary
                    V[ind] = overrelax[1]*V[ind] + overrelax[0]*(
                             a[ind][4]*V[ind-N[1]*N[2]] + a[ind][1]*V[ind+N[1]*N[2]] +
                             a[ind][5]*V[ind-N[2]] + a[ind][2]*V[ind+N[2]] +
                             a[ind][6]*V[ind-1] + a[ind][3]*V[ind+1] + rho[ind]*h2/EPS0)/a[ind][0]; //calculate new potential
                  } else { //no difference in permittivity, do not calculate new a values
                    V[ind] = overrelax[1]*V[ind] + overrelax[0]*(
                             V[ind-1*N[1]*N[2]] + V[ind+1*N[1]*N[2]] +
                             V[ind-1*N[2]] + V[ind+1*N[2]] +
                             V[ind-1] + V[ind+1] + rho[ind]*h2/EPS0/eps[ind]); //calculate new potential
                  }
                }
                if ( fabs((V[ind] - Vold)/ V[ind]) > currError){ //capture worst case error
                    currError = fabs((V[ind] - Vold)/V[ind]);
                }
              } else { //current cell IS electrode, set V[ind] = electrode voltage
                V[ind] = electrodemap[ind].second;
              }
            }
          }
        }
      }
      cycleCount++;
      if (cycleCount%50 == 0){
        std::cout << "On iteration " << cycleCount << " with " << currError*100 << "% error." << std::endl;
      }
  }while( currError > MAXERROR );
  std::cout << "Finished in " << cycleCount << " iterations." << std::endl;
  std::cout << "Time elapsed: " << float(clock()-begin_time)/CLOCKS_PER_SEC << " seconds" << std::endl;
  delete[] overrelax;
  delete[] isExterior;
  delete[] isBesideElec;
  delete[] isChangingeps;
  for ( int i = 1; i < N[0]-1; i++){ //for all x points except endpoints
    for ( int j = 1; j < N[1]-1; j++){
      for (int k = 1; k < N[2]-1; k++){
        ind = i*N[1]*N[2] + j*N[2] + k;
        delete[] a[ind];
      }
    }
  }
  delete[] a;
}

//Uses Successive Over Relaxation method to solve poisson's equation
void Solver::poisson3DSOR( void ){
  double Vold;
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
  static const std::vector<double> overrelax = {1.85/6, -0.85}; //overrelaxation parameter for SOR
  std::cout << "DOING SOR" << std::endl;
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
            Vold = V[ind];
            V[ind] = overrelax[1]*V[ind] + overrelax[0]*(
            V[ind-1*N[1]*N[2]] + V[ind+1*N[1]*N[2]] +
            V[ind-1*N[2]] + V[ind+1*N[2]] +
            V[ind-1] + V[ind+1] + rho[ind]*h2/EPS0); //calculate new potential
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

//Uses Iterative Jacobi method to solve poisson's equation
void Solver::poisson3DJacobi( void ){
  double* Vold = new double[sizeof(V)]; //needed to calculate error between new and old values
  double currError; //largest error on current loop
  unsigned int cycleCount = 0; //iteration count, for reporting
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
  delete[] Vold;
}

//Uses Iterative Gauss-Seidel method to solve poisson's equation
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
