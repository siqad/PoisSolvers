// @file:     phys_engine.h
// @author:   Samuel
// @created:  2017.08.23
// @editted:  2017.08.23 - Samuel
// @license:  GNU LGPL v3
//
// @desc:     Base class for physics engines

#ifndef _SIM_ANNEAL_PHYS_PHYS_ENGINE_H_
#define _SIM_ANNEAL_PHYS_PHYS_ENGINE_H_

#include "problem.h"
#include "electrodes.h"

#include <string>
#include <vector>
#include <boost/circular_buffer.hpp>



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

  namespace bpt = boost::property_tree;

  class PhysicsEngine
  {
  public:

    // constructor
    PhysicsEngine(const std::string &eng_nm, const std::string &i_path, const std::string &o_path);

    // destructor
    ~PhysicsEngine() {};

    // export results
    void writeResultsXml();

    // variables
    Problem problem;

    std::vector<Electrodes> elec_vec; // location of elecs
    double *arr; //will contain the potential

    std::vector<std::pair<float,float>> db_locs; // location of free dbs
    boost::circular_buffer<std::vector<int>> db_charges;

  private:
    std::string eng_name;
    std::string out_path;

  };

}

#endif
