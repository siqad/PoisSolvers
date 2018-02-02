#include "poissolver.h"

//Set the values in PhysConstants
namespace PhysConstants
{
   const double QE = 1.6021766208e-19; //elementary charge
   const double EPS0 = 8.85418782e-12; //permittivity of free space
   const double WF_GOLD = 5.1; //workfunction for gold in eV
   const double WF_COPPER = 4.7; //workfunction for copper in eV
   const double WF_ZINC = 4.3; //workfunction for zinc in eV
   const double WF_CESIUM = 2.1; //workfunction for cesium in eV
   const double WF_NICKEL = 5.01; //workfunction for nickel in eV
   const double CHI_SI = 4.05; //electron affinity for silicon in eV from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
   const double EPS_SI = 11.7; //relative permittivity of silicon
};


void PoisSolver::helloWorld(void){
  std::cout << "@@@@@@@@@@@@@HELLO WORLD@@@@@@@@@@@@" << std::endl;
}

void PoisSolver::worker(int step, std::vector<Electrodes> elec_vec){
  std::cout << "Modified by Nathan Chiu. Code is offered as is, with no warranty. See LICENCE.GPL for licence info." << std::endl;
  double cycleErr;
  int* indexErr = new int;

  // Calculate and save ds
  for (int i = 0; i<3; i++){
    SimParams::ds[i] = SimParams::Ls[i] / (double)SimParams::ns[i];
  }
  const int nsize = SimParams::ns[0]*SimParams::ns[1]*SimParams::ns[2]; //array size
  double *arr = new double[nsize];
  double *eps = new double[nsize]; //RELATIVE permittivity.
  double *chi = new double[nsize]; //electron affinity or WF
  double *RHS = new double[nsize]; // from which you can get a pointer to contiguous buffer, will contain rho.
  double *correction = new double[nsize]; // correction used to update RHS
  std::pair<int,double> *electrodemap = new std::pair<int,double>[nsize]; //stores electrode surface info and potentials.
  int cycleCount = 0;
  const std::clock_t begin_time = std::clock();
  init_eps(eps); // set permittivity
  init_rhs(chi, eps, RHS); // set rho and apply relative permittivity, also save electron affinity for bulk.
  create_electrode(RHS, electrodemap, chi, elec_vec);
  PoisFFT::Solver<3, double> S(SimParams::ns, SimParams::Ls, SimParams::BCs); //   create solver object, 3 dimensions, double precision

  std::cout << "Beginning solver" << std::endl;
  do{
    S.execute(arr, RHS); //run the solver, can be run many times for different right-hand side
    cycleErr = check_error(arr, correction, electrodemap, indexErr, eps);
    apply_correction(RHS, correction, electrodemap);
    std::cout << "On cycle " << cycleCount << " with error " << cycleErr << " at index " << *indexErr << ". " << arr[*indexErr] << " " << electrodemap[*indexErr].second << std::endl;
    cycleCount++;
  }while(cycleErr > SimParams::MAX_ERROR && cycleErr != 0);

  calc_charge(RHS, elec_vec);
  std::cout << "Finished on cycle " << cycleCount << " with error " << cycleErr << std::endl;

  save_file2D(arr, SimParams::OUTFILE); //solution is in arr
  save_file2D(RHS, SimParams::RHOFILE);
  save_file2D(correction, SimParams::CORRFILE);
  save_fileXML(arr, SimParams::RESXML, elec_vec);

  std::cout << "Time elapsed: " << float(clock()-begin_time)/CLOCKS_PER_SEC << " seconds" << std::endl;
  std::cout << "Ending, deleting variables" << std::endl;

  delete indexErr;
  delete[] eps;
  delete[] RHS;
  delete[] arr;
  delete[] electrodemap;
  delete[] chi;
  delete[] correction;
}

void PoisSolver::calc_charge(double* RHS , std::vector<Electrodes> elecs){
  //Want this to take in electrode location parameters, and spit out the charge on the conductor.
  double xmin, xmax, ymin, ymax, zmin, zmax, sum;
  for( int currElectrode = 0; currElectrode < elecs.size(); currElectrode++){
    xmin = elecs[currElectrode].x[0];
    xmax = elecs[currElectrode].x[1];
    ymin = elecs[currElectrode].y[0];
    ymax = elecs[currElectrode].y[1];
    zmin = elecs[currElectrode].z[0];
    zmax = elecs[currElectrode].z[1];
    sum = 0;
    for(int i = (int) SimParams::ns[0]*xmin/SimParams::Ls[0]; i <= (int) SimParams::ns[0]*xmax/SimParams::Ls[0]; i++){ //xmin first, then xmax
      for(int j = (int) SimParams::ns[1]*ymin/SimParams::Ls[1]; j <= (int) SimParams::ns[1]*ymax/SimParams::Ls[1]; j++){
        for(int k = (int) SimParams::ns[2]*zmin/SimParams::Ls[2]; k <= (int) SimParams::ns[2]*zmax/SimParams::Ls[2]; k++){
          sum += RHS[SimParams::IND(i,j,k)];
        }
      }
    }
    sum = sum*((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
    std::cout << "Calculated charge on electrode " << currElectrode << " = " << sum << std::endl;
  }
}

void PoisSolver::create_electrode(double* RHS, std::pair<int,double> *electrodemap, double* chi, std::vector<Electrodes> elecs){
  for(int i = 0; i < elecs.size(); i++){
    elecs[i].draw(SimParams::ns, SimParams::ds, SimParams::Ls, RHS, electrodemap, chi); //separately call draw for each electrode.
  }
}

void PoisSolver::apply_correction(double *RHS, double *correction, std::pair<int,double> *electrodemap){
  for(int i = 0; i < SimParams::ns[0]*SimParams::ns[1]*SimParams::ns[2]; i++){
    if(electrodemap[i].first == true){ //only correct error at electrode surfaces.
      RHS[i] -= correction[i];
    }
  }
}

double PoisSolver::check_error(double *arr, double *correction, std::pair<int,double> *electrodemap, int *indexErr, double *eps){
  double err = 0;
  double errOld;
  double correctionWeight;
  if(SimParams::BCs[0] == PoisFFT::PERIODIC){
    correctionWeight = 1e-7*SimParams::ns[0]*PhysConstants::EPS0/PhysConstants::QE/SimParams::Ls[0]/SimParams::Ls[0]; //Periodic
  } else if(SimParams::BCs[0] == PoisFFT::DIRICHLET || SimParams::BCs[0] == PoisFFT::NEUMANN){
    correctionWeight = 0.5e-7*SimParams::ns[0]*PhysConstants::EPS0/PhysConstants::QE/SimParams::Ls[0]/SimParams::Ls[0];  //Dirichlet & Neumann
  }
  for(int i = 0; i < SimParams::ns[0]*SimParams::ns[1]*SimParams::ns[2]; i++){
    if(electrodemap[i].first == true){ //only check error at electrodes.
      errOld = err;
      correction[i] = electrodemap[i].second-arr[i]; //intended potential - found potential. Looking at laplacian of potential doesn't allow snapping to electrode potentials.
      correction[i] *= correctionWeight;
      if(electrodemap[i].second != 0){
        err = std::max(err, fabs((arr[i] - electrodemap[i].second)/electrodemap[i].second)); //get largest error value.
      } else {
        err = std::max(err, arr[i]);
      }
      if( errOld != err ){
        *indexErr = i;
      }
    }
  }
  return err;
}

void PoisSolver::save_file2D(double* arr, char fname[]){
  std::string temp = fname;
  std::string finalpath = SimParams::resultpath + "/" + temp;
  std::ofstream outfile;
  outfile.open(finalpath.c_str(), std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping to " << finalpath << std::endl;
  const int k = SimParams::ns[2]/2;
  for (int i = 0; i < SimParams::ns[0]; i++){
    for (int j = 0; j < SimParams::ns[1]; j++){
      outfile << std::setprecision(5) << std::scientific << (i*SimParams::ds[0]-SimParams::xoffset)/SimParams::finalscale << " "
              << (j*SimParams::ds[1]-SimParams::yoffset)/SimParams::finalscale << " " << arr[SimParams::IND(i,j,k)] << std::endl;
    }
    outfile << std::endl;
  }
}


void PoisSolver::save_fileXML(double* arr, char fname[], std::vector<Electrodes> elec_vec){

  // define major XML nodes
  boost::property_tree::ptree tree;
  boost::property_tree::ptree node_root;       // <sim_out>
  boost::property_tree::ptree node_eng_info;   // <eng_info>
  boost::property_tree::ptree node_sim_params; // <sim_params>
  boost::property_tree::ptree node_electrode;  // <electrode>
  boost::property_tree::ptree node_potential_map;  // <potential>
  std::string temp = fname;
  std::string finalpath = SimParams::resultpath + "/" + temp;

  std::cout << "Write results to XML..." << std::endl;
  // NOTE in the future, there's probably a range of stuff that can be exported.
  // for now, only export charge config

  // eng_info
  node_eng_info.put("engine", "PoisSolver");
  node_eng_info.put("version", "TBD"); // TODO real version

  // sim_params
  // TODO

  // electrode
  for (auto elec : elec_vec) {
    boost::property_tree::ptree node_dim;
    node_dim.put("<xmlattr>.x1", (std::to_string((elec.x[0]-SimParams::xoffset)/SimParams::finalscale/SimParams::Ls[0]).c_str()));
    node_dim.put("<xmlattr>.y1", (std::to_string((elec.y[0]-SimParams::yoffset)/SimParams::finalscale/SimParams::Ls[1]).c_str()));
    node_dim.put("<xmlattr>.x2", (std::to_string((elec.x[1]-SimParams::xoffset)/SimParams::finalscale/SimParams::Ls[0]).c_str()));
    node_dim.put("<xmlattr>.y2", (std::to_string((elec.y[1]-SimParams::yoffset)/SimParams::finalscale/SimParams::Ls[1]).c_str()));
    node_electrode.add_child("dim", node_dim);
    boost::property_tree::ptree node_pot;
    node_pot.put("", std::to_string(elec.potential).c_str());
    node_electrode.add_child("potential", node_pot);
  }

  //potential_map
  const int k = SimParams::ns[2]/2;
  for (int i = 0; i < SimParams::ns[0]; i++){
    for (int j = 0; j < SimParams::ns[1]; j++){
      //create each entry
      boost::property_tree::ptree node_potential_val;
      node_potential_val.put("<xmlattr>.x", (std::to_string((i*SimParams::ds[0]-SimParams::xoffset)/SimParams::finalscale/SimParams::Ls[0])).c_str());
      node_potential_val.put("<xmlattr>.y", (std::to_string((j*SimParams::ds[1]-SimParams::yoffset)/SimParams::finalscale/SimParams::Ls[1])).c_str());
      node_potential_val.put("<xmlattr>.val", std::to_string(arr[SimParams::IND(i,j,k)]).c_str());
      node_potential_map.add_child("potential_val", node_potential_val);
    }
  }
  node_root.add_child("eng_info", node_eng_info);
  node_root.add_child("electrode", node_electrode);
  node_root.add_child("potential_map", node_potential_map);
  tree.add_child("sim_out", node_root);

  // write to file
  boost::property_tree::write_xml(finalpath, tree, std::locale(), boost::property_tree::xml_writer_make_settings<std::string>(' ',4));

  std::cout << "Write to XML complete." << std::endl;
}


void PoisSolver::init_eps(double* eps){
  int i,j,k;
  std::cout << "Initialising eps" << std::endl;
  for (i=0;i<SimParams::ns[0];i++){
    double x = SimParams::ds[0]*(i+0.5);
    for (j=0;j<SimParams::ns[1];j++){
      double y = SimParams::ds[1]*(j+0.5);
      for (k=0;k<SimParams::ns[2];k++){
        double z = SimParams::ds[2]*(k+0.5);
        if (y < SimParams::Ls[1]/2){
          // a[IND(i,j,k)] = PhysConstants::EPS_SI; //Si relative permittivity
          eps[SimParams::IND(i,j,k)] = 1;  //Free space
        } else{
          eps[SimParams::IND(i,j,k)] = 1;  //Free space
        }
      }
    }
  }
  std::cout << "Finished eps initialisation" << std::endl;
}

void PoisSolver::init_rhs(double* chi, double* eps, double* rhs){
  int i,j,k;
  std::cout << "Initialising RHS" << std::endl;
  for (i=0;i<SimParams::ns[0];i++){
    double x = SimParams::ds[0]*(i+0.5);
    for (j=0;j<SimParams::ns[1];j++){
      double y = SimParams::ds[1]*(j+0.5);
      for (k=0;k<SimParams::ns[2];k++){
        double z = SimParams::ds[2]*(k+0.5);
        // a[IND(i,j,k)] = 1e16*PhysConstants::QE/PhysConstants::EPS0/eps[IND(i,j,k)]; //in m^-3, scale by permittivity
        rhs[SimParams::IND(i,j,k)] = 0; //in m^-3, scale by permittivity
        chi[SimParams::IND(i,j,k)] = PhysConstants::CHI_SI;
      }
    }
  }
  std::cout << "Finished RHS initialisation" << std::endl;
}
