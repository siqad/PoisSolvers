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
#include "electrodes.hpp"
#include <vector>

const double pi = 3.14159265358979323846;

#define EPS0 8.85418782e-12
#define Q_E 1.6e-19
#define MAXERROR 1e-1
#define IND(i,j,k) (i)*(ns[1]*ns[2])+(j)*(ns[2])+k
#define FILENAME ((char*)"outfile.txt")
#define RHOFILE ((char*) "outrho.txt")
#define EPSFILE ((char*) "outeps.txt")
#define CORRECTIONFILE ((char*) "outcorr.txt")
#define WF_GOLD 5.1 //workfunction for gold in eV
#define WF_COPPER 4.7 //workfunction for copper in eV
#define WF_ZINC 4.3 //workfunction for zinc in eV
#define WF_CESIUM 2.1 //workfunction for cesium in eV
#define WF_NICKEL 5.01
#define CHI_SI 4.05//electron affinity for silicon in eV from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
#define EPS_SI 11.7

void save_file2D(const int[], const double[], const double[], double*, char[]);
void check_solution(const int[], const double[], const double[], double*);
void init_eps(const int[], const double[], const double[], double*);
void init_rhs(const int[], const double[], const double[], double*, double*, double*);
double check_error(const int[], const double[], double*, double*, std::pair<int,double>*, int*, double*, double*, const int[]);
void apply_correction(const int[], const double[], double*, double*, std::pair<int,double>*);
void create_electrode(const int[], const double[], const double[], double*, std::pair<int,double>*, double*, const int, Electrodes[]);
void worker(int, std::vector<Electrodes>);
void calc_charge(const int[], const double[], double*, const int, Electrodes[]);
void temp(int*, std::vector<Electrodes>*);

int main(void){
  int count = 0;
  std::vector<Electrodes> elec_vec;
  temp(&count, &elec_vec);
  std::cout << "VECTOR SIZE " << elec_vec.size() << std::endl;
  for(int i = 0; i < 1; i++){
    worker(i, elec_vec);
  }
  return 0;
}

void temp(int *count, std::vector<Electrodes> *elecs){

  int pix_x1, pix_x2, pix_y1, pix_y2;
  double potential;
  boost::property_tree::ptree tree; // Create empty property tree object
  boost::property_tree::read_xml("cooldbdesign.xml", tree); // Parse the XML into the property tree.
  BOOST_FOREACH(boost::property_tree::ptree::value_type &node, tree.get_child("dbdesigner.design")) {
    boost::property_tree::ptree subtree = node.second; //get subtree with layer items at the top
    if( node.first == "layer"){ //go one level below layers.
      BOOST_FOREACH( boost::property_tree::ptree::value_type const&v, subtree.get_child( "" ) ) {
        boost::property_tree::ptree subtree2 = v.second; //get subtree with layer item params at the top
        if (v.first == "electrode"){
          BOOST_FOREACH( boost::property_tree::ptree::value_type const&v2, subtree2.get_child( "" ) ) {
            std::string label = v2.first; //get the name of each param
            if(label == "dim"){ //Read the 4 numbers
              *count = *count + 1;
              pix_x1 = v2.second.get<int>("<xmlattr>.x1", -1); //get the 2D corners in pixel distances
              pix_x2 = v2.second.get<int>("<xmlattr>.x2", -1);
              pix_y1 = v2.second.get<int>("<xmlattr>.y1", -1);
              pix_y2 = v2.second.get<int>("<xmlattr>.y2", -1);
              std::cout << pix_x1 << " " << pix_y1 << ", " << pix_x2 << " " << pix_y2 << std::endl;
              // elecs->push_back(Electrodes(pix_x1*1e-8, pix_x2*1e-8, pix_y1*1e-8, pix_y2*1e-8, 0.3e-6, 0.7e-6, potential, WF_GOLD));
            }else if(label == "potential") {
              potential = subtree2.get<double>(label);
              std::cout << label << ":  " << potential << std::endl;
            }else if(label == "electrode_type"){ //electrode_type is the last one
              elecs->push_back(Electrodes(pix_x1*1e-8, pix_x2*1e-8, pix_y1*1e-8, pix_y2*1e-8, 0.3e-6, 0.7e-6, potential, WF_GOLD));
            }else if(label != "<xmlattr>"){ //unexpected extras
              std::string value = subtree2.get<std::string>(label);
              std::cout << label << ":  " << value << std::endl;
            }
          }
        }
      }
    }
  }
}

void worker(int step, std::vector<Electrodes> elec_vec){
  std::cout << "Modified by Nathan Chiu. Code is offered as is, with no warranty. See LICENCE.GPL for licence info." << std::endl;
  const double Ls[3] = {1.0e-6, 1.0e-6, 1.0e-6}; //x, y, z domain dimensions in MICROMETRE
  const int ns[3] = {100, 100, 100}; //x, y, z gridpoint numbers
  double ds[3];  // distances between gridpoints
  double cycleErr;
  int* indexErr = new int;
  const int BCs[6] = {PoisFFT::NEUMANN, PoisFFT::NEUMANN,  //boundary conditions
                      PoisFFT::NEUMANN, PoisFFT::NEUMANN,
                      PoisFFT::NEUMANN, PoisFFT::NEUMANN};
  int i;
  for (i = 0; i<3; i++){ // set the grid, depends on the boundary conditions
    ds[i] = Ls[i] / ns[i];
  }
  // const int numElectrodes = 4;
  // const int numElectrodes = 2;
  const int numElectrodes = elec_vec.size();
  Electrodes elecs[numElectrodes];
  for (int i = 0; i < numElectrodes; i++){
    elecs[i] = elec_vec[i];
    std::cout << elecs[i].x[0] << std::endl;
    std::cout << elecs[i].x[1] << std::endl;
    std::cout << elecs[i].y[0] << std::endl;
    std::cout << elecs[i].y[1] << std::endl;
    std::cout << elecs[i].z[0] << std::endl;
    std::cout << elecs[i].z[1] << std::endl;
    std::cout << elecs[i].potential << std::endl;
    std::cout << elecs[i].WF << std::endl;
  }
  // Electrodes elecs[numElectrodes] = {
  //   Electrodes(0.1e-6, 0.37e-6, 0.04e-6, 0.38e-6, 0.3e-6, 0.7e-6, 6, WF_GOLD),
  //   Electrodes(0.28e-6, 0.74e-6, 0.55e-6, 0.84e-6, 0.3e-6, 0.7e-6, 5, WF_GOLD),
  //   // Electrodes(0.6e-6, 0.8e-6, 0.2e-6, 0.4e-6, 0.3e-6, 0.7e-6, 1, WF_GOLD),
  //   // Electrodes(0.6e-6, 0.8e-6, 0.6e-6, 0.8e-6, 0.3e-6, 0.7e-6, 1, WF_GOLD)
  // };
  elecs[step].potential = 2.0;
  double *arr = new double[ns[0]*ns[1]*ns[2]];
  double *arrOld = new double[ns[0]*ns[1]*ns[2]];
  double *eps = new double[ns[0]*ns[1]*ns[2]]; //RELATIVE permittivity.
  double *chi = new double[ns[0]*ns[1]*ns[2]]; //electron affinity or WF
  double *RHS = new double[ns[0]*ns[1]*ns[2]]; // from which you can get a pointer to contiguous buffer, will contain rho.
  double *correction = new double[ns[0]*ns[1]*ns[2]]; // correction used to update RHS
  std::pair<int,double> *electrodemap = new std::pair<int,double>[ns[0]*ns[1]*ns[2]]; //stores electrode surface info and potentials.
  int cycleCount = 0;
  const std::clock_t begin_time = std::clock();
  init_eps(ns, ds, Ls, eps); // set permittivity
  init_rhs(ns, ds, Ls, chi, eps, RHS); // set rho and apply relative permittivity, also save electron affinity for bulk.
  create_electrode(ns, ds, Ls, RHS, electrodemap, chi, numElectrodes, elecs);
  PoisFFT::Solver<3, double> S(ns, Ls, BCs); //   create solver object, 3 dimensions, double precision
  std::cout << "Beginning solver" << std::endl;
  ProfilerStart("./profileresult.out"); //using google performance tools
  do{
    S.execute(arr, RHS); //run the solver, can be run many times for different right-hand side
    cycleErr = check_error(ns, Ls, arr, correction, electrodemap, indexErr, arrOld, eps, BCs);
    apply_correction(ns, ds, RHS, correction, electrodemap);
    std::cout << "On cycle " << cycleCount << " with error " << cycleErr << " at index " << *indexErr << ". " << arr[*indexErr] << " " << electrodemap[*indexErr].second << std::endl;
    cycleCount++;
  }while(cycleErr > MAXERROR && cycleErr != 0);
  calc_charge(ns, Ls, RHS, numElectrodes, elecs);
  ProfilerStop();
  std::cout << "Finished on cycle " << cycleCount << " with error " << cycleErr << std::endl;
  save_file2D(ns, ds, Ls, RHS, RHOFILE);
  save_file2D(ns, ds, Ls, arr, FILENAME); //solution is in arr
  save_file2D(ns, ds, Ls, correction, CORRECTIONFILE);
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

void calc_charge( const int ns[3], const double Ls[3], double* RHS , const int numElectrodes, Electrodes elecs[]){
  //Want this to take in electrode location parameters, and spit out the charge on the conductor.
  double xmin, xmax, ymin, ymax, zmin, zmax, sum;
  for( int currElectrode = 0; currElectrode < numElectrodes; currElectrode++){
    xmin = elecs[currElectrode].x[0];
    xmax = elecs[currElectrode].x[1];
    ymin = elecs[currElectrode].y[0];
    ymax = elecs[currElectrode].y[1];
    zmin = elecs[currElectrode].z[0];
    zmax = elecs[currElectrode].z[1];
    sum = 0;
    for(int i = (int) ns[0]*xmin/Ls[0]; i <= (int) ns[0]*xmax/Ls[0]; i++){ //xmin first, then xmax
      for(int j = (int) ns[1]*ymin/Ls[1]; j <= (int) ns[1]*ymax/Ls[1]; j++){
        for(int k = (int) ns[2]*zmin/Ls[2]; k <= (int) ns[2]*zmax/Ls[2]; k++){
          sum += RHS[IND(i,j,k)];
        }
      }
    }
    sum = sum*((xmax-xmin)*(ymax-ymin)*(zmax-zmin));
    std::cout << "Calculated charge on electrode " << currElectrode << " = " << sum << std::endl;
  }
}

void create_electrode(const int ns[3], const double ds[3], const double Ls[3], double* RHS, std::pair<int,double> *electrodemap, double* chi, const int numElectrodes, Electrodes elecs[]){
  for(int i = 0; i < numElectrodes; i++){
    elecs[i].draw(ns, ds, Ls, RHS, electrodemap, chi); //separately call draw for each electrode.
  }
}

void apply_correction(const int ns[3], const double ds[3], double *RHS, double *correction, std::pair<int,double> *electrodemap){
  for(int i = 0; i < ns[0]*ns[1]*ns[2]; i++){
    if(electrodemap[i].first == true){ //only correct error at electrode surfaces.
      RHS[i] -= correction[i];
    }
  }
}

double check_error(const int ns[3], const double Ls[3], double *arr, double *correction, std::pair<int,double> *electrodemap, int *indexErr,
                   double *arrOld, double *eps, const int BCs[6]){
  double err = 0;
  double errOld;
  double correctionWeight;
  if(BCs[0] == PoisFFT::PERIODIC){
    correctionWeight = 1e-7*ns[0]*EPS0/Q_E/Ls[0]/Ls[0]; //Periodic
  } else if(BCs[0] == PoisFFT::DIRICHLET || BCs[0] == PoisFFT::NEUMANN){
    correctionWeight = 0.5e-7*ns[0]*EPS0/Q_E/Ls[0]/Ls[0];  //Dirichlet & Neumann
  }
  // for(int i = 0; i < ns[0]*ns[1]*ns[2]; i++){
  for(int i = 0; i < ns[0]*ns[1]*ns[2]; i++){
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
    // arrOld[i] = arr[i];
  }
  return err;
}

void save_file2D(const int ns[3], const double ds[3], const double Ls[3], double* arr, char fname[]){
  std::ofstream outfile;
  outfile.open(fname, std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping to " << fname << std::endl;
  const int k = ns[2]/2;
  for (int i = 0; i < ns[0]; i++){
    for (int j = 0; j < ns[1]; j++){
        outfile << std::setprecision(5) << std::scientific << i*ds[0] << " " << j*ds[1] << " " << arr[IND(i,j,k)] << std::endl;
    }
    outfile << std::endl;
  }
}

void init_eps(const int ns[3], const double ds[3], const double Ls[3], double* a){
  int i,j,k;
  std::cout << "Initialising eps" << std::endl;
  for (i=0;i<ns[0];i++){
    double x = ds[0]*(i+0.5);
    for (j=0;j<ns[1];j++){
      double y = ds[1]*(j+0.5);
      for (k=0;k<ns[2];k++){
        double z = ds[2]*(k+0.5);
        if (y < Ls[1]/2){
          // a[IND(i,j,k)] = EPS_SI; //Si relative permittivity
          a[IND(i,j,k)] = 1;  //Free space
        } else{
          // a[IND(i,j,k)] = EPS_SI; //Si relative permittivity
          a[IND(i,j,k)] = 1;  //Free space
        }
      }
    }
  }
  std::cout << "Finished rho initialisation" << std::endl;
}

void init_rhs(const int ns[3], const double ds[3], const double Ls[3], double* chi, double* eps, double* a){
  int i,j,k;
  std::cout << "Initialising RHS" << std::endl;
  for (i=0;i<ns[0];i++){
    double x = ds[0]*(i+0.5);
    for (j=0;j<ns[1];j++){
      double y = ds[1]*(j+0.5);
      for (k=0;k<ns[2];k++){
        double z = ds[2]*(k+0.5);
        // a[IND(i,j,k)] = 1e16*Q_E/EPS0/eps[IND(i,j,k)]; //in m^-3, scale by permittivity
        a[IND(i,j,k)] = 0; //in m^-3, scale by permittivity
        chi[IND(i,j,k)] = CHI_SI;
      }
    }
  }
  std::cout << "Finished RHS initialisation" << std::endl;
}
