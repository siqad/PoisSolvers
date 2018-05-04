#include "poissolver.h"
#include "python2.7/Python.h"
#include <math.h>
#include <stdio.h>
#include <Python.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

namespace SimParams
{
  // scaling and offset values, populated in setBuffer
  double finalscale, xoffset, yoffset, zoffset;
  std::string resultpath;

  //stuff used during simulation
  double Ls[3] = {1e-6, 1e-6, 1e-6}; // simulation length in x, y, z
  int ns[3] = {50, 50, 50}; // resolution in x, y, z
  double ds[3]; // simulation length per resolution step, CALCULATED
  // int BCs[6]  = {PoisFFT::NEUMANN, PoisFFT::NEUMANN,
  //                PoisFFT::NEUMANN, PoisFFT::NEUMANN,
  //                PoisFFT::NEUMANN, PoisFFT::NEUMANN}; // boundary conditions for left, right, top, bottom, front, back.
  int BCs[6]  = {0, 0,
                 0, 0,
                 0, 0}; // boundary conditions for left, right, top, bottom, front, back.
  double MAX_ERROR = 5e-2;
  int IND(int i, int j, int k){ return (i)*(ns[1]*ns[2]) + (j)*(ns[2]) + k; };

  //stuff used post-simulation
  char* RESXML = (char*) "sim_result.xml";
};


using namespace phys;


int main(){
  // std::cout << "Number of command line arguments: " << argc << std::endl;
  // std::vector<Electrodes> elec_vec;
  // std::string in_path;
  // std::string out_path;
  
  std::cout << "HELLO" << std::endl;
  std::cout << "Calling Python with new protocol" << std::endl;
  // Py_Initialize();
  // std::string script_path = "./test.py";
  // PyObject *obj = Py_BuildValue("s", script_path.c_str());
  // FILE *file = _Py_fopen( "test.py", "r+" );
  // // FILE* script_file = _Py_fopen_obj(obj, "r+");
  // //FILE* script_file = fopen(scriptPath().c_str(), "r");
  // int argc = 0;
  // wchar_t * argv[1];
  // 
  // argv[0] = Py_DecodeLocale(script_path.c_str(), NULL);
  // 
  // PyRun_SimpleFile(script_file, script_path.c_str());
  // Py_Finalize();
  // std::string script_path = "./test.py";
  int argc = 3;
  // wchar_t * argv[2];
  char* argv[3];
  // wchar_t* a;
  // stringToWChar(a, a1);
  // wchar_t* wpath;
  // stringToWChar(wpath, script_path);
  // argv[0] = Py_DecodeLocale(script_path.c_str(), NULL);
  // argv[0] = Py_DecodeLocale(scriptPath().c_str(), NULL);
  // argv[1] = Py_DecodeLocale("-i", NULL);
  // argv[2] = Py_DecodeLocale(script_problem_path.c_str(), NULL);
  // argv[3] = Py_DecodeLocale("-o", NULL);
  // argv[4] = Py_DecodeLocale(script_result_path.c_str(), NULL);

  std::string script_name = "test.py";
  argv[0] = &script_name[0u];
  std::string in_path = "hello1" ;
  argv[1] = &in_path[0u];
  std::string out_path = "hello2" ;
  argv[2] = &out_path[0u];
  
  // std::cout << a1 << std::endl;
  // argv[1] = a;
  std::string script_path = "./test.py";
  FILE* PythonScriptFile = fopen(script_path.c_str(), "r");
  if(PythonScriptFile)
  {
    Py_SetProgramName(argv[0]);
    std::cout <<  Py_GetProgramName() << std::endl;
    Py_Initialize();
    PySys_SetArgv(argc, argv);
    PyRun_SimpleFile(PythonScriptFile, script_path.c_str());
    fclose(PythonScriptFile);
  }

  Py_Finalize();




  return 0;
}

// int stringToWChar(wchar_t* &wc, const std::string &s)
// {
//     std::wstring wsTmp(s.begin(), s.end());
//     wc = wsTmp.c_str();
//     return 0;
// }
