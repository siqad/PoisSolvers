#include "poissolver.h"

namespace SimParams
{
  // scaling and offset values, populated in set_buffer
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


using namespace phys;


int main(int argc,char* argv[]){
  std::cout << "Number of command line arguments: " << argc << std::endl;
  std::vector<Electrodes> elec_vec;
  std::string in_path;
  std::string out_path;

  if(argc < 3){
    std::cout << "Insufficient inputs provided, program terminating." << std::endl;
  }else{ //need at LEAST binary call, input path, output path.
    in_path = argv[1];
    out_path = argv[2];
    // std::cout << SimParams::resultpath << std::endl;


    std::cout << std::endl << "*** Constructing Problem ***" << std::endl;
    PoisSolver ps(in_path, out_path);
    ps.initVars();

    std::cout << std::endl << "*** Run Simulation ***" << std::endl;
    if(!ps.runSim()) {
      std::cout << "Simulation failed, aborting" << std::endl;
      return 0;
    }
    //writing to XML taken care of inside poissolver

  }
  return 0;
}
