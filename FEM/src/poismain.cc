#include "poissolver.h"

using namespace phys;

int main(int argc,char* argv[]){

  std::string in_path, out_path;
  if (argc > 2) {
    // Assuming first 2 command line args are sim_problem location and sim_result location
    in_path = argv[1];
    out_path = argv[2];
    std::cout << "Input: " << in_path << std::endl;
    std::cout << "Output: " << out_path << std::endl;
    std::cout << std::endl << "*** Constructing Problem ***" << std::endl;
    PoisSolver ps(in_path, out_path);
    ps.runSolver();
  } else {
    std::cout << "Call to program must be in the form /path/to/poissolver /path/to/input /path/to/output" << std::endl;
  }
  return 0;
}
