#include "poissolver.h"
#include "Python.h"
#include <stdio.h>

using namespace phys;

PoisSolver::PoisSolver(const std::string& i_path, const std::string& o_path)
{
  in_path = i_path;
  out_path = o_path;
  std::cout << "PoisSolver initialized." << std::endl;
}

void PoisSolver::runSolver(void)
{
  std::cout << "Running Solver" << std::endl;
  int argc = 4;
  wchar_t * argv[4];
  std::string script_name = "Poisson_3D.py";
  argv[0] = Py_DecodeLocale(&script_name[0u], NULL);
  argv[1] = Py_DecodeLocale(&in_path[0u], NULL);
  argv[2] = Py_DecodeLocale(&out_path[0u], NULL);
  std::string working_dir = "./build/debug/src/phys/poissolver/";
  argv[3] = Py_DecodeLocale(&working_dir[0u], NULL);

  std::string script_path = "./build/debug/src/phys/poissolver/Poisson_3D.py";
  FILE* PythonScriptFile = fopen(script_path.c_str(), "r+");
  if(PythonScriptFile)
  {

    std::cout << "*** Calling Python with new protocol ***" << std::endl;
    Py_SetProgramName(argv[0]);
    Py_Initialize();
    PySys_SetArgv(argc, argv);
    PyRun_SimpleFile(PythonScriptFile, script_path.c_str());
    fclose(PythonScriptFile);
    Py_Finalize();
  }
}
