#include "poissolver.h"

using namespace phys;

int main(int argc,char* argv[]){
  
  std::string in_path, out_path;
  if (argc > 2) {
    // Assuming first 2 command line args are sim_problem location and sim_result location
    in_path = argv[1];
    out_path = argv[2];
  } else {
    in_path = "../sim_problem.xml" ;
    out_path = "../sim_result.xml" ;
  }
  
  std::cout << "Input: " << in_path << std::endl;
  std::cout << "Output: " << out_path << std::endl;
  std::cout << std::endl << "*** Constructing Problem ***" << std::endl;
  PoisSolver ps(in_path, out_path);
  ps.initVars();
  ps.runSolver();
  return 0;
}