#include <utility>
#include "electrodes.h"
#include <iostream>

Electrodes::Electrodes( void ){}


Electrodes::~Electrodes( void ){}


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

int Electrodes::IND(int i, int j, int k, const int ns[]){
  return (i)*(ns[1]*ns[2]) + (j)*(ns[2]) + k;
}

void Electrodes::draw(const int ns[3], const double Ls[3], double* RHS, std::pair<int,double> *electrodemap, double* chi){
  int i, j, k; //draw the electrode into an electrode map
  // set RHS inside electrodes to be 0, E field should be 0 inside metal.
  for(i = (int) ns[0]*x[0]/Ls[0]; i < (int) ns[0]*x[1]/Ls[0]; i++){
    for(j = (int) ns[1]*y[0]/Ls[1]; j < (int) ns[1]*y[1]/Ls[1]; j++){
      for(k = (int) ns[2]*z[0]/Ls[2]; k < (int) ns[2]*z[1]/Ls[2]; k++){
        RHS[IND(i,j,k,ns)] = 0;
      }
    }
  }

  // set metal and immediate surrounding semiconductor electrodemap values
  // to be true, to indicate a fixed potential. Then, set the potential.
  // Repeat for all dimensions.
  for( int iter = 0; iter < 2; iter++){
    i = (int) ns[0]*x[iter]/Ls[0]; //xmin first, then xmax
    for(j = (int) ns[1]*y[0]/Ls[1]; j <= (int) ns[1]*y[1]/Ls[1]; j++){
      for(k = (int) ns[2]*z[0]/Ls[2]; k <= (int) ns[2]*z[1]/Ls[2]; k++){
        electrodemap[IND(i,j,k,ns)].first = true;
        electrodemap[IND(i,j,k,ns)].second = potential;
        if((x[0]!=0 && iter==0) || (x[1]!=ns[0]-1 && iter==1)){
          electrodemap[IND(i-1+2*iter,j,k,ns)].first = true;
          electrodemap[IND(i-1+2*iter,j,k,ns)].second = potential - (WF-chi[IND(i-1+2*iter,j,k,ns)]);
        }
      }
    }
    j = (int) ns[1]*y[iter]/Ls[1]; //ymin first, then ymax
    for(i = (int) ns[0]*x[0]/Ls[0]; i <= (int) ns[0]*x[1]/Ls[0]; i++){
      for(k = (int) ns[2]*z[0]/Ls[2]; k <= (int) ns[2]*z[1]/Ls[2]; k++){
        electrodemap[IND(i,j,k,ns)].first = true; //set true for electrode surface.
        electrodemap[IND(i,j,k,ns)].second = potential; //set electrode potential
        if((y[0]!=0 && iter==0) || (y[1]!=ns[1]-1 && iter==1)){ //set potential of adjacent silicon with WF.
          electrodemap[IND(i,j-1+2*iter,k,ns)].first = true;
          electrodemap[IND(i,j-1+2*iter,k,ns)].second = potential - (WF-chi[IND(i,j-1+2*iter,k,ns)]);
        }
      }
    }
    k = (int) ns[2]*z[iter]/Ls[2]; //zmin
    for(i = (int) ns[0]*x[0]/Ls[0]; i <= (int) ns[0]*x[1]/Ls[0]; i++){
      for(j = (int) ns[1]*y[0]/Ls[1]; j <= (int) ns[1]*y[1]/Ls[1]; j++){
        electrodemap[IND(i,j,k,ns)].first = true; //set true for electrode surface.
        electrodemap[IND(i,j,k,ns)].second = potential; //set electrode potential
        if((z[0]!=0 && iter==0) || (z[1]!=ns[2]-1 && iter==1)){ //set potential of adjacent silicon with WF.
          electrodemap[IND(i,j,k-1+2*iter,ns)].first = true;
          electrodemap[IND(i,j,k-1+2*iter,ns)].second = potential - (WF-chi[IND(i,j,k-1+2*iter,ns)]);
        }
      }
    }
  }
}
