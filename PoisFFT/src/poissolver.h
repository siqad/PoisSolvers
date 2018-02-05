#ifndef POISSOLVER_H
#define POISSOLVER_H

#include <cmath>
#include <iostream>
#include "poisfft.h"
#include <stdlib.h>
#include <iomanip>
#include <fstream>
#include <utility>
#include <algorithm>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include "electrodes.h"
#include <vector>
#include <string>
#include "phys_engine.h"

// Physical constants used for calculation
namespace PhysConstants
{
   extern const double QE; //elementary charge
   extern const double EPS0; //permittivity of free space
   extern const double WF_GOLD; //workfunction for gold in eV
   extern const double WF_COPPER; //workfunction for copper in eV
   extern const double WF_ZINC; //workfunction for zinc in eV
   extern const double WF_CESIUM; //workfunction for cesium in eV
   extern const double WF_NICKEL; //workfunction for nickel in eV
   extern const double CHI_SI; //electron affinity for silicon in eV from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
   extern const double EPS_SI; //relative permittivity of silicon
};

// Simulation parameters used for simulation control
namespace SimParams
{
  // scaling and offset values
  extern double finalscale, xoffset, yoffset;
  extern std::string resultpath;

  //stuff used during simulation
  extern double Ls[3]; // simulation length in x, y, z
  extern int ns[3]; // resolution in x, y, z
  extern double ds[3]; // simulation length per resolution step, CALCULATED
  extern int BCs[6]; // boundary conditions for left, right, top, bottom, front, back.
  extern double MAX_ERROR;
  extern int IND(int i, int j, int k);

  //stuff used post-simulation
  extern char* OUTFILE;
  extern char* RHOFILE;
  extern char* EPSFILE;
  extern char* CORRFILE;
  extern char* RESXML;
};


namespace phys{
  class PoisSolver: public PhysicsEngine
  {
    public:
        PoisSolver(const std::string& i_path, const std::string& o_path);
        ~PoisSolver(){} //destructor1
        void helloWorld(void);
        void worker(int step, std::vector<Electrodes> elec_vec);
        void init_eps(double* eps);
        void init_rhs(double* chi, double* eps, double* rhs);
        void save_file2D(double* arr, char fname[]);
        void save_fileXML(double* arr, char fname[], std::vector<Electrodes> elec_vec);
        double check_error(double *arr, double *correction, std::pair<int,double> *electrodemap, int *indexErr, double *eps);
        void apply_correction(double *RHS, double *correction, std::pair<int,double> *electrodemap);
        void create_electrode(double* RHS, std::pair<int,double> *electrodemap, double* chi, std::vector<Electrodes> elecs);
        void calc_charge(double* RHS , std::vector<Electrodes> elecs);
        bool runSim();
        void rescale();
        std::vector<Electrodes> set_buffer(std::vector<Electrodes> elec_vec);
  };
}
// functions used by main
void parse_tree(std::vector<Electrodes> *elecs, std::string path);

#endif //POISSOLVER_H
