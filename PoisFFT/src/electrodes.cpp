#include <utility>
#include "electrodes.h"


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

  for(int i = (int) x[0]/Ls[0]*ns[0]; i < (int) x[1]/Ls[0]*ns[0]; i++){ //set RHS 0 inside electrodes
    for(int j = (int) y[0]/Ls[1]*ns[1]; j < (int) y[1]/Ls[1]*ns[1]; j++){
      for(int k = (int) z[0]/Ls[2]*ns[2]; k < (int) z[1]/Ls[2]*ns[2]; k++){
        RHS[i*ns[1]*ns[2] + j*ns[2] + k] = 0;
      }
    }
  }
  for( int iter = 0; iter < 2; iter++){
    i = (int) x[iter]/Ls[0]*ns[0]; //xmin first, then xmax
    for(int j = (int) y[0]/Ls[1]*ns[1]; j < (int) y[1]/Ls[1]*ns[1]; j++){
      for(int k = (int) z[0]/Ls[2]*ns[2]; k < (int) z[1]/Ls[2]*ns[2]; k++){
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].first = true; //set true for electrode surface.
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].second = potential; //set electrode potential
      }
    }
  }
  for( int iter = 0; iter < 2; iter++){
    j = (int) y[iter]/Ls[1]*ns[1]; //ymin first, then ymax
    for(int i = (int) x[0]/Ls[0]*ns[0]; i < (int) x[1]/Ls[0]*ns[0]; i++){
      for(int k = (int) z[0]/Ls[2]*ns[2]; k < (int) z[1]/Ls[2]*ns[2]; k++){
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].first = true; //set true for electrode surface.
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].second = potential; //set electrode potential
      }
    }
  }
  for( int iter = 0; iter < 2; iter++){
    k = (int) z[iter]/Ls[1]*ns[1]; //zmin
    for(int i = (int) x[0]/Ls[0]*ns[0]; i < (int) x[1]/Ls[0]*ns[0]; i++){
      for(int j = (int) y[0]/Ls[1]*ns[1]; j < (int) y[1]/Ls[1]*ns[1]; j++){
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].first = true; //set true for electrode surface.
        electrodemap[i*ns[1]*ns[2] + j*ns[2] + k].second = potential; //set electrode potential
      }
    }
  }
}
