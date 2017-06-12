#include "electrodes.h"
#include "solver.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <ctime>
#include <fstream>
#include <cmath>

Electrodes::Electrodes( void ) {
  dims.resize(3);
  centre.resize(3);
}
Electrodes::~Electrodes( void ) {
  delete[] pdims;
  delete[] pcentre;
}

//initialize sizes
void Electrodes::init_elec( void ){
 //initialized to zero.

}

//draw the electrode into solver's electrode map
void Electrodes::draw( Solver * s ){
  //dimensions must be odd, in order to have equal lengths on left and right sides
  //for now, will round up to next even number
  for( int i = pcentre[0] - pdims[0]/2; i < pcentre[0] + pdims[0]/2; i++){
    for( int j = pcentre[1] - pdims[1]/2; j < pcentre[1] + pdims[1]/2; j++){
      for( int k = pcentre[2] - pdims[2]/2; k < pcentre[2] + pdims[2]/2; k++){
        s->electrodemap[i*s->N[1]*s->N[2] + j*s->N[2] + k].first = workfunction; //set electrode work function
        s->electrodemap[i*s->N[1]*s->N[2] + j*s->N[2] + k].second = potential; //set electrode potential
      }
    }
  }
}

void Electrodes::convert( void ){
  pdims = &dims[0];
  pcentre = &centre[0];
}
