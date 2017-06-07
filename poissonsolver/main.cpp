 #include <iostream>
#include "./src/solver.h"
#include "./src/electrodes.h"

int main() {
    //Parameters
    int Nx = 100;
    int Ny = 100;
    int Nz = 100;
    double Lx = 1.0;
    double Ly = 1.0;
    double Lz = 1.0;
    int numelectrodes = 4;
    Solver sol;
    Solver * s = &sol;
    std::vector<Electrodes> elec(numelectrodes);
    s->N = s->set_val(Nx, Ny, Nz);
    s->set_val(s->L, Lx, Ly, Lz);
//REQUIRE CONSISTENT L/N IN ALL X, Y, Z
    s->h2 = Lx*Lx/Nx/Nx;
    s->h = Lx/Nx;
    s->set_val(s->boundarytype, NEUMANN);
    //initialize
    for( int index = 0; index < numelectrodes; index++){
      elec[index].init_elec(); //init all electrodes
    }
    //4 total electrodes, arranged in a square
    elec[0].centre = {Nx/4, Ny/4, Nz/2};
    elec[0].dims = {Nx/8, Ny/8, Nz/8};
    elec[0].potential = 0;
    elec[0].workfunction = WF_COPPER;
    elec[1].centre = {3*Nx/4, 3*Ny/4, Nz/2};
    elec[1].dims = {Nx/8, Ny/8, Nz/8};
    elec[1].potential = 1e-10;
    elec[1].workfunction = WF_COPPER;
    elec[2].centre = {3*Nx/4, Ny/4, Nz/2};
    elec[2].dims = {Nx/8, Ny/8, Nz/8};
    elec[2].potential = 2e-10;
    elec[2].workfunction = WF_GOLD;
    elec[3].centre = {Nx/4, 3*Ny/4, Nz/2};
    elec[3].dims = {Nx/8, Ny/8, Nz/8};
    elec[3].potential = 4e-10;
    elec[3].workfunction = WF_GOLD;
    s->rho = s->init_val( 0, s->rho);
    s->V = s->init_val( 0, s->V ); //Need to initiate V before setting Boundary conds.
    s->electrodemap = s->init_val( 0, s-> electrodemap );
    for( int index = 0; index < numelectrodes; index++){
      elec[index].draw(s); //draw all electrodes into electrodemap
    }

    s->init_rho( );
    s->init_eps( ); //will replace with reading from file eventually
    s->set_BCs(0, 0, 0, 0, 0, 0); //boundary conditions
    //reset solution vector and call for SOR_GEN
    s->set_val( s->solvemethod, SOR_GEN);
    s->solve();
    s->write_2D(s->V, FILENAMESOR_GEN);
    s->write_2D(s->rho, FILENAMERHO);

//need to delete all the arrays created with new. (init_val() and set_val())
std::cout << "DELETING VARIABLES" << std::endl;
    s->del(s->V);
    s->del(s->N);
    s->del(s->rho);
    s->del(s->electrodemap);
    return 0;
}
