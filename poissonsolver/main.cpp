 #include <iostream>
#include "./src/solver.h"
#include "./src/electrodes.h"

int main() {
    //Parameters
    int Nx = 100;
    int Ny = 100;
    int Nz = 100;
    double Lx = 1;
    double Ly = 2;
    double Lz = 3;
    int numelectrodes = 2;
    Solver sol;
    Solver * s = &sol;
    std::vector<Electrodes> elec(numelectrodes);
    s->N = s->set_val(Nx, Ny, Nz);
    s->set_val(s->L, Lx, Ly, Lz);
    s->set_val(s->h2, LATTICECONSTANT*LATTICECONSTANT);
    s->set_val(s->h, LATTICECONSTANT);
    s->set_val(s->boundarytype, NEUMANN);
    //initialize
    for( int index = 0; index < numelectrodes; index++){
      elec[index].init_elec(); //init all electrodes
    }
    elec[0].centre = {Nx/4, Ny/4, Nz/2};
    elec[0].dims = {Nx/8, Ny/8, Nz/8};
    elec[1].centre = {3*Nx/4, 3*Ny/4, Nz/2};
    elec[1].dims = {Nx/8, Ny/8, Nz/8};
    s->rho = s->init_val( 0, s->rho);
    s->V = s->init_val( 0, s->V ); //Need to initiate V before setting Boundary conds.
    s->electrodemap = s->init_val( 0, s-> electrodemap );
    for( int index = 0; index < numelectrodes; index++){
      elec[index].draw(s); //draw all electrodes into electrodemap
    }

    s->init_rho( );
    s->init_eps( ); //will replace with reading from file eventually
    s->set_BCs(0, 0, 0, 0, 0, 0); //Dirichlet boundary conditions
    //reset solution vector and call for SOR_GEN
    s->set_val( s->solvemethod, SOR_GEN);
    s->solve();
    s->write_2D(s->V, FILENAMESOR_GEN);
    s->write_2D(s->rho, FILENAMERHO);
    s->del_V();
    delete[] s->N;
    return 0;
}
