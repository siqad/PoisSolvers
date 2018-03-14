#include <iostream>
#include "./src/solver.h"
#include "./src/electrodes.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <string>
#include <cmath>

void show_array(const boost::numeric::ublas::unbounded_array<double>& a)
{
    for(size_t i=0; i<a.size(); ++i)
      std::cout << a[i] << ' ';
    std::cout << '\n';
}

void workingcall( double t ){
//Parameters
  int Nx = 10;
  int Ny = 10;
  int Nz = 10;
  double Lx = 1.0;
  double Ly = 1.0;
  double Lz = 1.0;
  int numelectrodes = 6;
  int torig = t;
  Solver * s = new Solver();
  std::vector<Electrodes> elec(numelectrodes);
  std::cout << "torig =" << torig << std::endl;
  t = t/64;
  std::cout << "timestep:" << t << std::endl;
      s->set_val(Nx, Ny, Nz, s->N);
      s->set_val(Lx, Ly, Lz, s->L);
      s->h2 = Lx*Lx/Nx/Nx;
      s->h = Lx/Nx;
      s->boundarytype = NEUMANN;
  //can allocate space for rho, eps, and V now that number of lattice points is known.
      s->rho = s->init_val( 0, s->rho);
      s->V = s->init_val( 0, s->V );
      s->eps = s->init_val( 0, s->eps );
  //REQUIRE CONSISTENT L/N IN ALL X, Y, Z
      s->init_rho( );
      s->initEPS( );
  //initialize electrodes
      for( int index = 0; index < numelectrodes; index++){
        elec[index].init_elec(); //init all electrodes
      }
  //6 total electrodes, arranged in a square
      elec[0].centre = {Nx*2/100, Ny/2, Nz/4};
      elec[0].dims = {Nx*4/100, Ny/2, Nz/11};
      // elec[0].potential = 10;
      elec[0].potential = 10*sin(-(2*PI*t));
      elec[0].workfunction = WF_GOLD;
      elec[1].centre = {Nx*25/100, Ny/2, Nz/4};
      elec[1].dims = {Nx/10, Ny/2, Nz/11};
      // elec[1].potential = -10;
      elec[1].potential = 10*sin(-(2*PI*t+PI/2));
      elec[1].workfunction = WF_GOLD;
      elec[2].centre = {Nx*50/100, Ny/2, Nz/4};
      elec[2].dims = {Nx/10, Ny/2, Nz/11};
      // elec[2].potential = 10;
      elec[2].potential = 10*sin(-(2*PI*t+PI));
      elec[2].workfunction = WF_GOLD;
      elec[3].centre = {Nx*75/100, Ny/2, Nz/4};
      elec[3].dims = {Nx/10, Ny/2, Nz/11};
      // elec[3].potential = -10;
      elec[3].potential = 10*sin(-(2*PI*t+3*PI/2));
      elec[3].workfunction = WF_GOLD;
      elec[4].centre = {Nx*98/100, Ny/2, Nz/4};
      elec[4].dims = {Nx*4/100, Ny/2, Nz/11};
      // elec[4].potential = 10;
      elec[4].potential = 10*sin(-(2*PI*t));
      elec[4].workfunction = WF_GOLD;
      elec[5].centre = {Nx/2, Ny/2, (int) 9.5*Nz/10};
      elec[5].dims = {Nx, Ny, Nz/10};
      // elec[5].potential = 0;
      elec[5].potential = 0;
      elec[5].workfunction = CHI_SI;

      s->electrodemap = s->init_val( 0, s-> electrodemap );
      for( int index = 0; index < numelectrodes; index++){
        elec[index].convert(); //convert vectors into pointers
        elec[index].draw(s); //draw all electrodes into electrodemap
      }
  //eps and rho are already pointers. To use existing rho and eps, initialize and set with:
      s->set_BCs(0, 0, 0, 0, 0, 0); //boundary conditions
  //reset solution vector and call for SOR_GEN
      s->solvemethod = SOR_BLAS;
      s->solve();

      s->write_2D(s->V, FILENAMEBLAS + std::to_string( (int) torig ) + ".txt");
      s->write_2D(s->rho, FILENAMERHO);
      s->write_2D(s->eps, FILENAMEEPS);
  delete s;
//destructors take care of deleting.
}

int main() {
    // for (int i = 0; i < 64; i++){
    //   workingcall( (double) i );
    // }
    workingcall(0);
    return 0;
}
