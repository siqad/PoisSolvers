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
    Solver sol;
    Electrodes elec0;
    Electrodes elec1;

    //initialize
    Solver * s = &sol;
    Electrodes * pelec = &elec0;
    s->set_val(s->N, Nx, Ny, Nz);
    s->set_val(s->L, Lx, Ly, Lz);
    s->set_val(s->h2, LATTICECONSTANT*LATTICECONSTANT);
    s->init_rho( );
    s->init_eps( ); //will replace with reading from file eventually
    s->init_val( s->V, 0); //Need to initiate V before setting Boundary conds.
    s->set_BCs(0, 0, 0, 0, 0, 0); //Dirichlet boundary conditions
    pelec->init_elec(s);
    pelec->centre = {Nx/4, Ny/4, Nz/2};
    pelec->dims = {Nx/8, Ny/8, Nz/8};
    pelec->draw(s);
    pelec->centre = {3*Nx/4, 3*Ny/4, Nz/2};
    pelec->dims = {Nx/8, Ny/8, Nz/8};
    pelec->draw(s);

    //reset solution vector and call for SOR_GEN
    s->init_val( s->V, 0);
    s->set_val( s->solvemethod, SOR_GEN);
    s->solve();
    s->write_2D(s->V, FILENAMESOR_GEN);
    s->write_2D(s->rho, FILENAMERHO);
    return 0;
}
