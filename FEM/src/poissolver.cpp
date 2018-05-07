#include "poissolver.h"
// #include <Python.h>
#include "python2.7/Python.h"
#include <stdio.h>

using namespace phys;

PoisSolver::PoisSolver(const std::string& i_path, const std::string& o_path)
{
  phys_con = new PhysicsConnector(std::string("PoisSolver"), i_path, o_path);
  initSolver();
}

void PoisSolver::initSolver(void)
{
  std::cout << "PoisSolver instantiated." << std::endl;
  phys_con->setRequiredSimParam("bcs");
  phys_con->setRequiredSimParam("resolution");
  // phys_con->setRequiredSimParam("length");
  phys_con->setRequiredSimParam("max_abs_error");
  phys_con->setRequiredSimParam("max_rel_error");
  phys_con->setRequiredSimParam("mode");
  phys_con->setRequiredSimParam("high_pot");
  phys_con->setRequiredSimParam("low_pot");
  phys_con->setRequiredSimParam("steps");
  phys_con->setExpectElectrode(true);
  phys_con->setExpectDB(true);
  phys_con->readProblem();
  for (auto& iter : phys_con->getRequiredSimParam()) {
    if(!phys_con->parameterExists(iter)){
      std::cout << "Parameter " << iter << " not found." << std::endl;
    }
  }
}


void PoisSolver::initVars(void)
{
  std::cout << "PoisSolver::initVars" << std::endl;
  //pulling the required sim params from PhysicsConnector
  bc = phys_con->parameterExists("bcs") ?
                  phys_con->getParameter("bcs") : "robin";
  resolution = phys_con->parameterExists("resolution") ?
                  std::stof(phys_con->getParameter("resolution")) : 1;
  // length = phys_con->parameterExists("length") ?
  //                 std::stod(phys_con->getParameter("length")) : 1e-6;
  max_error = phys_con->parameterExists("max_error") ?
                  std::stod(phys_con->getParameter("max_error")) : 5e-2;
  mode = phys_con->parameterExists("mode") ?
                  phys_con->getParameter("mode") : "Standard";
  high_pot = phys_con->parameterExists("high_pot") ?
                  std::stod(phys_con->getParameter("high_pot")) : 0;
  low_pot = phys_con->parameterExists("low_pot") ?
                  std::stod(phys_con->getParameter("low_pot")) : 0;
  steps = phys_con->parameterExists("steps") ?
                std::stod(phys_con->getParameter("steps")) : 20;

  std::map<std::string, std::string> metal_props = phys_con->getProperty("Metal");

  //metal layer properties
  z_offset = std::stod((std::string)(metal_props["zoffset"]));
  z_thickness = std::stod((std::string)(metal_props["zheight"]));
  
  std::cout << "Parameters initialized" << std::endl;
}


void PoisSolver::runSolver(void)
{
  std::cout << "Running Solver" << std::endl;
  int argc = 4;
  char* argv[4];
  std::string script_name = "Poisson_3D.py";
  // std::string script_name = "./build/debug/src/phys/poissolver";
  argv[0] = &script_name[0u];
  std::string in_path = phys_con->getInputPath();
  argv[1] = &in_path[0u];
  std::string out_path = phys_con->getOutputPath();
  argv[2] = &out_path[0u];
  std::string working_dir = "./build/debug/src/phys/poissolver/";
  argv[3] = &working_dir[0u];
  
  float wf_gold = 5.1;
  phys_con->initCollections();
  
  std::cout << "Saving electrodes" << std::endl;
  for(auto elec : *(phys_con->elec_col)) {
    elec_vec.push_back(Electrodes(
      elec->x1, elec->x2,
      elec->y1, elec->y2,
      z_offset, z_offset + z_thickness,
      elec->potential, wf_gold)); //elec_vec is part of phys_engine
  }
  
  
  // std::string script_path = "Poisson_3D.py";
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
  loadPotentialTxt(e_pot_2D);
  exportData();
  // std::cout << e_pot_2D.size() << std::endl;
}

void PoisSolver::loadPotentialTxt(std::vector<std::vector<float>> &elec_data)
{
  boost::filesystem::path p(phys_con->getInputPath());
  boost::filesystem::path dir = p.parent_path();
  dir /= "arrays_tmp.txt";
  std::ifstream inFile;
  inFile.open(dir.string());
  float x, y, z;
  // std::vector<std::vector<float>> elec_data;
  elec_data.clear();
  while(inFile >> x >> y >> z){
    std::vector<float> v{x,y,z};
    elec_data.push_back(v);
  }
}

void PoisSolver::exportData(void){
  //create the vector of [x, y, val] that will be sent to setElecPotentialData
  //number of data points is n*n
  std::vector<std::vector<std::string>> elec_pot_data(e_pot_2D.size());
  for (unsigned int i = 0; i < e_pot_2D.size(); i++){
    elec_pot_data[i].resize(3);
    elec_pot_data[i][0] = std::to_string(e_pot_2D[i][0]);
    elec_pot_data[i][1] = std::to_string(e_pot_2D[i][1]);
    elec_pot_data[i][2] = std::to_string(e_pot_2D[i][2]);
  }
  phys_con->setExportElecPotential(true);
  phys_con->setElecPotentialData(elec_pot_data);

  std::vector<std::vector<std::string>> electrode_data(elec_vec.size());
  double x1, x2, y1, y2, pot;
  for (unsigned int i = 0; i < elec_vec.size(); i++){
    electrode_data[i].resize(5);
    x1 = elec_vec[i].x[0];
    x2 = elec_vec[i].x[1];
    y1 = elec_vec[i].y[0];
    y2 = elec_vec[i].y[1];
    pot = elec_vec[i].potential;
    electrode_data[i][0] = std::to_string(x1);
    electrode_data[i][1] = std::to_string(y1);
    electrode_data[i][2] = std::to_string(x2);
    electrode_data[i][3] = std::to_string(y2);
    electrode_data[i][4] = std::to_string(pot);
  }

  phys_con->setExportElectrode(true);
  phys_con->setElectrodeData(electrode_data);

  std::cout << std::endl << "*** Write Result to Output ***" << std::endl;
  phys_con->writeResultsXml();

}
