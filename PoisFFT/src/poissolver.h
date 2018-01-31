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
#include <gperftools/profiler.h>
#include <ctime>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include "electrodes.h"
#include <vector>
#include <string>


// Physical constants used for calculation
namespace PhysConstants
{
   const double QE = 1.6021766208e-19; //elementary charge
  //  const double PI = 3.14159265358979323846;
   const double EPS0 = 8.85418782e-12; //permittivity of free space
   const double WF_GOLD = 5.1; //workfunction for gold in eV
   const double WF_COPPER = 4.7; //workfunction for copper in eV
   const double WF_ZINC = 4.3; //workfunction for zinc in eV
   const double WF_CESIUM = 2.1; //workfunction for cesium in eV
   const double WF_NICKEL = 5.01; //workfunction for nickel in eV
   const double CHI_SI = 4.05; //electron affinity for silicon in eV from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
   const double EPS_SI = 11.7; //relative permittivity of silicon
};

// Simulation parameters used for simulation control
// ds needs to be reinitialised in source, since it depends on Ls and ns.
namespace SimParams
{
  // scaling and offset values
  double finalscale, xoffset, yoffset;
  std::string resultpath;

  //stuff used during simulation
  double Ls[3] = {1e-6, 1e-6, 1e-6}; // simulation length in x, y, z
  int ns[3] = {50, 50, 50}; // resolution in x, y, z
  double ds[3]; // simulation length per resolution step, CALCULATED
  int BCs[6]  = {PoisFFT::NEUMANN, PoisFFT::NEUMANN,
                 PoisFFT::NEUMANN, PoisFFT::NEUMANN,
                 PoisFFT::NEUMANN, PoisFFT::NEUMANN}; // boundary conditions for left, right, top, bottom, front, back.
  double MAX_ERROR = 5e-2;
  int IND(int i, int j, int k){ return (i)*(ns[1]*ns[2]) + (j)*(ns[2]) + k; };

  //stuff used post-simulation
  char* OUTFILE = (char*) "outfile.txt";
  char* RHOFILE = (char*) "rhofile.txt";
  char* EPSFILE = (char*) "epsfile.txt";
  char* CORRFILE = (char*) "corrfile.txt";
  char* RESXML = (char*) "sim_result.xml";
};

class Poissolver
{
  // public:
  //   Electrodes(); //default constructor
  //   Electrodes(double, double, double, double, double, double, double, double); //parametrized constructor
  //   ~Electrodes(); //destructor
  //   double x[2]; //xmin (x[0]) and xmax (x[1])
  //   double y[2]; //ymin and ymax
  //   double z[2]; //zmin and zmax
  //   double potential;   //pointer after conversion of vector
  //   double WF;
  //   void draw(const int[3], const double[3], const double[3], double*, std::pair<int,double>*, double*);
};

void save_file2D(double* arr, char fname[]);
void save_fileXML(double* arr, char fname[], std::vector<Electrodes> elec_vec);
void init_eps(double* eps);
void init_rhs(double* chi, double* eps, double* rhs);
double check_error(double *arr, double *correction, std::pair<int,double> *electrodemap, int *indexErr, double *eps);
void apply_correction(double *RHS, double *correction, std::pair<int,double> *electrodemap);
void create_electrode(double* RHS, std::pair<int,double> *electrodemap, double* chi, std::vector<Electrodes> elecs);
void worker(int step, std::vector<Electrodes> elec_vec);
void calc_charge(double* RHS , std::vector<Electrodes> elecs);
void parse_tree(std::vector<Electrodes> *elecs, std::string path);
std::vector<Electrodes> set_buffer(std::vector<Electrodes> elec_vec);

#endif //POISSOLVER_H
