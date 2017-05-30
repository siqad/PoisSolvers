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

    s->set_val(s->N, Nx, Ny, Nz);
    s->set_val(s->L, Lx, Ly, Lz);
    s->set_val(s->h2, LATTICECONSTANT*LATTICECONSTANT);
    //initialize
    for( int index = 0; index < numelectrodes; index++){
      elec[index].init_elec(); //init all electrodes
    }
    elec[0].centre = {Nx/4, Ny/4, Nz/2};
    elec[0].dims = {Nx/8, Ny/8, Nz/8};
    elec[1].centre = {3*Nx/4, 3*Ny/4, Nz/2};
    elec[1].dims = {Nx/8, Ny/8, Nz/8};
    s->electrodemap.resize(Nx*Ny*Nz); //size electrodemap
    for( int index = 0; index < numelectrodes; index++){
      elec[index].draw(s);
    }

    s->init_rho( );
    s->init_eps( ); //will replace with reading from file eventually
    s->init_val( s->V, 0); //Need to initiate V before setting Boundary conds.
    s->set_BCs(0, 0, 0, 0, 0, 0); //Dirichlet boundary conditions


    //reset solution vector and call for SOR_GEN
    s->init_val( s->V, 0);
    s->set_val( s->solvemethod, SOR_GEN);
    s->solve();
    s->write_2D(s->V, FILENAMESOR_GEN);
    s->write_2D(s->rho, FILENAMERHO);
    return 0;
}
