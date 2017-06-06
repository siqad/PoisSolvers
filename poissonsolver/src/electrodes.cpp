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

//initialize sizes
void Electrodes::init_elec( void ){
 //initialized to zero.
  dims.resize(3);
  centre.resize(3);
}

//draw the electrode into solver's electrode map
void Electrodes::draw( Solver * s ){
  //dimensions must be odd, in order to have equal lengths on left and right sides
  //for now, will round up to next even number

  for( int i = centre[0] - dims[0]/2; i < centre[0] + dims[0]/2; i++){
    for( int j = centre[1] - dims[1]/2; j < centre[1] + dims[1]/2; j++){
      for( int k = centre[2] - dims[2]/2; k < centre[2] + dims[2]/2; k++){
        s->electrodemap[i*s->N[1]*s->N[2] + j*s->N[2] + k].first = workfunction; //set electrode work function
        s->electrodemap[i*s->N[1]*s->N[2] + j*s->N[2] + k].second = potential; //set electrode potential
      }
    }
  }
}
