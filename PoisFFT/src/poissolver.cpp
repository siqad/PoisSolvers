#include <cmath>
#include <iostream>
#include "poisfft.h"
#include "poissolver.h"
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

// Physical constants are placed inside PhysConstants namespace.

#define SIMGRID 50
#define SIMLENGTH 1.0e-6
#define MAXERROR 5e-2
#define IND(i,j,k) (i)*(ns[1]*ns[2])+(j)*(ns[2])+k
#define FILENAME ((char*)"outfile.txt")
#define RHOFILE ((char*) "outrho.txt")
#define EPSFILE ((char*) "outeps.txt")
#define CORRECTIONFILE ((char*) "outcorr.txt")


void save_file2D(const int[], const double[], const double[], double*, char[], std::string, double, double, double);
void save_fileXML(const int[], const double[], const double[], double*, char[], std::string, double, double, double, std::vector<Electrodes>);

void check_solution(const int[], const double[], const double[], double*);
void init_eps(const int[], const double[], const double[], double*);
void init_rhs(const int[], const double[], const double[], double*, double*, double*);
double check_error(const int[], const double[], double*, double*, std::pair<int,double>*, int*, double*, double*, const int[]);
void apply_correction(const int[], const double[], double*, double*, std::pair<int,double>*);
void create_electrode(const int[], const double[], const double[], double*, std::pair<int,double>*, double*, const int, Electrodes[]);

void worker(int, std::vector<Electrodes>, std::string, double, double, double);
void calc_charge(const int[], const double[], double*, const int, Electrodes[]);
void parse_tree(std::vector<Electrodes>*, std::string);
std::vector<Electrodes> set_buffer(std::vector<Electrodes>, double*, double*, double*);

int main(int argc,char* argv[]){

  std::cout << "Number of command line arguments: " << argc << std::endl;
  std::vector<Electrodes> elec_vec;
  std::string arg1;
  std::string arg2;
  std::string outpath;
  double finalscale = 0;
  double xoffset = 0;
  double yoffset = 0;

  for(int i = 0; i < argc; i++){
    std::cout << argv[i] << std::endl;
  }

  if(argc == 1){
    std::cout << "No path was passed to the solver program. Program terminating." << std::endl;
  }else if(argc > 2){ //need at LEAST binary call, input path, output path.
    arg1 = argv[1];
    arg2 = argv[2];
    if(arg1.find(".xml") != std::string::npos){ //argv[1] is an xml path, assume it is the INPUT file.
      std::cout << "Input path detected." << std::endl;
      parse_tree(&elec_vec, argv[1]);
      elec_vec = set_buffer(elec_vec, &finalscale, &xoffset, &yoffset);
      std::cout << finalscale << " " << xoffset << " " << yoffset << std::endl;
      std::cout << "VECTOR SIZE " << elec_vec.size() << std::endl;
      if(arg2.find(".xml") != std::string::npos){ //argv[2] is an xml path, assume it is the OUTPUT file.
        std::cout << "Output path detected." << std::endl;
        std::size_t found = arg2.find_last_of("/\\");
        outpath = arg2.substr(0, found);
        std::cout << outpath << std::endl;
      }else{
        std::cout << "Path not detected. Result XML path needs to be provided as second argument to binary. Terminating." << std::endl;
      }
      for(int i = 0; i < 1; i++){
        worker(i, elec_vec, outpath, finalscale, xoffset, yoffset); //where the magic happens
      }

    }else{
      std::cout << "Path not detected. Problem XML path needs to be provided as first argument to binary. Terminating." << std::endl;
    }
  }

  return 0;
}

std::vector<Electrodes> set_buffer(std::vector<Electrodes> elec_vec, double* finalscale, double* xoffset, double* yoffset) {
  //want to scale electrodes down to fit the simulation space.
  //First, find the (min, max) (x, y) values.
  //Then, extend or contract so that 10% of the sim space is left on each edge as buffer space.
  double xmin = 100;
  double ymin = 100;
  double xmax = 0;
  double ymax = 0;
  double xlength;
  double ylength;
  double xscale = 1;
  double yscale = 1;
  // double finalscale;
  // double xoffset = 0;
  // double yoffset = 0;
  for(int i = 0; i < elec_vec.size(); i++){
    xmin = std::min(xmin, elec_vec[i].x[0]);
    ymin = std::min(ymin, elec_vec[i].y[0]);
    xmax = std::max(xmax, elec_vec[i].x[1]);
    ymax = std::max(ymax, elec_vec[i].y[1]);
  }

  std::cout << "xmin: " << xmin << std::endl;
  std::cout << "xmax: " << xmax << std::endl;
  std::cout << "ymin: " << ymin << std::endl;
  std::cout << "ymax: " << ymax << std::endl;

  //x-scaling to keep 10% buffer on each horizontal side.
  if(xmin < 0.1*SIMLENGTH || xmax > 0.9*SIMLENGTH){
    xlength = xmax - xmin;
    xscale = 0.8*SIMLENGTH/xlength;
    std::cout << "xlength: " << xlength << std::endl;
  }
  if(ymin < 0.1*SIMLENGTH || ymax > 0.9*SIMLENGTH){
    ylength = ymax - ymin;
    yscale = 0.8*SIMLENGTH/ylength;
    std::cout << "ylength: " << ylength << std::endl;
  }
  // std::cout << "xscale: " << xscale << std::endl;
  // std::cout << "yscale: " << yscale << std::endl;
  //scale all elements by lowest scaling factor.
  *finalscale = std::min(xscale, yscale);
  std::cout << "Final scaling factor is: " << *finalscale << std::endl;
  for(int i = 0; i < elec_vec.size(); i++){
    std::cout << elec_vec[i].x[0] << std::endl;
    std::cout << elec_vec[i].x[1] << std::endl;
    std::cout << elec_vec[i].y[0] << std::endl;
    std::cout << elec_vec[i].y[1] << std::endl;

    elec_vec[i].x[0] *= *(finalscale);
    elec_vec[i].x[1] *= *(finalscale);
    elec_vec[i].y[0] *= *(finalscale);
    elec_vec[i].y[1] *= *(finalscale);

    std::cout << elec_vec[i].x[0] << std::endl;
    std::cout << elec_vec[i].x[1] << std::endl;
    std::cout << elec_vec[i].y[0] << std::endl;
    std::cout << elec_vec[i].y[1] << std::endl;

  }
  //now sample is sure to fit within simulation boundaries, with space for buffer.
  //translate the violating part to the buffer boundary, once for x and once for y.
  xmin = 100;
  ymin = 100;
  xmax = 0;
  ymax = 0;
  //find how far outside the boundary the shapes still sit.
  for(int i = 0; i < elec_vec.size(); i++){
    xmin = std::min(xmin, elec_vec[i].x[0]);
    ymin = std::min(ymin, elec_vec[i].y[0]);
    xmax = std::max(xmax, elec_vec[i].x[1]);
    ymax = std::max(ymax, elec_vec[i].y[1]);
  }
  std::cout << "xmin: " << xmin << std::endl;
  std::cout << "xmax: " << xmax << std::endl;
  std::cout << "ymin: " << ymin << std::endl;
  std::cout << "ymax: " << ymax << std::endl;

  if(xmin < 0.1*SIMLENGTH){  //too far to the left, want positive offset to bring it right
    //find the offset
    *xoffset = 0.1*SIMLENGTH - xmin;
  }else if(xmax > 0.9*SIMLENGTH){ //too far right, want negative offset to bring it left.
    // *xoffset = xmax - 0.9*SIMLENGTH;
    *xoffset = 0.9*SIMLENGTH - xmax;
  }
  if(ymin < 0.1*SIMLENGTH){ //too far up
    //find the offset in y
    *yoffset = 0.1*SIMLENGTH - ymin;
  }else if(ymax > 0.9*SIMLENGTH){ //too far down
    // *yoffset = ymax - 0.9*SIMLENGTH;
    *yoffset = 0.9*SIMLENGTH - ymax;
  }

  std::cout << "X offset is: " << *xoffset << std::endl;
  std::cout << "Y offset is: " << *yoffset << std::endl;
  //fix the offsets
  for(int i = 0; i < elec_vec.size(); i++){ //move all points based on offset.
    elec_vec[i].x[0] += *(xoffset);
    elec_vec[i].x[1] += *(xoffset);
    elec_vec[i].y[0] += *(yoffset);
    elec_vec[i].y[1] += *(yoffset);
    // std::cout << elec_vec[i].x[0] << " " << elec_vec[i].x[1] << ", " << elec_vec[i].y[0] << " " << elec_vec[i].y[1] << std::endl;
  }
  return elec_vec;

}

void parse_tree(std::vector<Electrodes> *elecs, std::string path){

  int pix_x1, pix_x2, pix_y1, pix_y2;
  double potential;
  std::cout << "PARSING NOW, PATH NAME IS " << path << std::endl;
  boost::property_tree::ptree tree; // Create empty property tree object
  // boost::property_tree::read_xml("cooldbdesign.xml", tree); // Parse the XML into the property tree.
  boost::property_tree::read_xml(path, tree); // Parse the XML into the property tree.
  BOOST_FOREACH(boost::property_tree::ptree::value_type &node, tree.get_child("dbdesigner.design")) {
    boost::property_tree::ptree subtree = node.second; //get subtree with layer items at the top
    if( node.first == "layer"){ //go one level below layers.
      std::string type = node.second.get<std::string>("<xmlattr>.type");
      if (type == "Electrode"){ //make sure that the layer type is Electrode
        BOOST_FOREACH( boost::property_tree::ptree::value_type const&v, subtree.get_child( "" ) ) {
          boost::property_tree::ptree subtree2 = v.second; //get subtree with layer item params at the top
          if (v.first == "electrode"){ //for each electrode, read into memory
            BOOST_FOREACH( boost::property_tree::ptree::value_type const&v2, subtree2.get_child( "" ) ) {
              std::string label = v2.first; //get the name of each param
              if(label == "dim"){ //Read the 4 numbers
                pix_x1 = v2.second.get<int>("<xmlattr>.x1", -1); //get the 2D corners in pixel distances
                pix_x2 = v2.second.get<int>("<xmlattr>.x2", -1);
                pix_y1 = v2.second.get<int>("<xmlattr>.y1", -1);
                pix_y2 = v2.second.get<int>("<xmlattr>.y2", -1);
                std::cout << pix_x1 << " " << pix_y1 << ", " << pix_x2 << " " << pix_y2 << std::endl;
                // elecs->push_back(Electrodes(pix_x1*1e-8, pix_x2*1e-8, pix_y1*1e-8, pix_y2*1e-8, 0.3e-6, 0.7e-6, potential, PhysConstants::WF_GOLD));
              }else if(label == "potential") {
                potential = subtree2.get<double>(label);
                std::cout << label << ":  " << potential << std::endl;
              }else if(label == "electrode_type"){ //electrode_type is the last one
                elecs->push_back(Electrodes(pix_x1*SIMLENGTH, pix_x2*SIMLENGTH, pix_y1*SIMLENGTH, pix_y2*SIMLENGTH, 0.3e-6, 0.7e-6, potential, PhysConstants::WF_GOLD));
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
  std::cout << "Successfully read " << path << std::endl;
}


void worker(int step, std::vector<Electrodes> elec_vec, std::string resultpath, double finalscale, double xoffset, double yoffset){
  std::cout << "Modified by Nathan Chiu. Code is offered as is, with no warranty. See LICENCE.GPL for licence info." << std::endl;


  const double Ls[3] = {SIMLENGTH, SIMLENGTH, SIMLENGTH}; //x, y, z domain dimensions in MICROMETRE
  const int ns[3] = {SIMGRID, SIMGRID, SIMGRID}; //x, y, z gridpoint numbers
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
  do{
    S.execute(arr, RHS); //run the solver, can be run many times for different right-hand side
    cycleErr = check_error(ns, Ls, arr, correction, electrodemap, indexErr, arrOld, eps, BCs);
    apply_correction(ns, ds, RHS, correction, electrodemap);
    std::cout << "On cycle " << cycleCount << " with error " << cycleErr << " at index " << *indexErr << ". " << arr[*indexErr] << " " << electrodemap[*indexErr].second << std::endl;
    cycleCount++;
  }while(cycleErr > MAXERROR && cycleErr != 0);


  calc_charge(ns, Ls, RHS, numElectrodes, elecs);
  std::cout << "Finished on cycle " << cycleCount << " with error " << cycleErr << std::endl;


  save_file2D(ns, ds, Ls, RHS, RHOFILE, resultpath, finalscale, xoffset, yoffset);
  save_file2D(ns, ds, Ls, arr, FILENAME, resultpath, finalscale, xoffset, yoffset); //solution is in arr
  save_file2D(ns, ds, Ls, correction, CORRECTIONFILE, resultpath, finalscale, xoffset, yoffset);
  save_fileXML(ns, ds, Ls, arr, "sim_result.xml", resultpath, finalscale, xoffset, yoffset, elec_vec);


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
    correctionWeight = 1e-7*ns[0]*PhysConstants::EPS0/PhysConstants::QE/Ls[0]/Ls[0]; //Periodic
  } else if(BCs[0] == PoisFFT::DIRICHLET || BCs[0] == PoisFFT::NEUMANN){
    correctionWeight = 0.5e-7*ns[0]*PhysConstants::EPS0/PhysConstants::QE/Ls[0]/Ls[0];  //Dirichlet & Neumann
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


void save_file2D(const int ns[3], const double ds[3], const double Ls[3], double* arr, char fname[], std::string pathname, double finalscale, double xoffset, double yoffset){
  // std::string finalpath = fname;
  std::string temp = fname;
  std::string finalpath = pathname+"/"+temp;
  std::ofstream outfile;
  // std::cout << pathname << std::endl;
  // std::cout << finalpath << std::endl;
  // std::cout << temp << std::endl;
  outfile.open(finalpath.c_str(), std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping to " << finalpath << std::endl;
  const int k = ns[2]/2;
  for (int i = 0; i < ns[0]; i++){
    for (int j = 0; j < ns[1]; j++){
      outfile << std::setprecision(5) << std::scientific << (i*ds[0]-xoffset)/finalscale << " "
              << (j*ds[1]-yoffset)/finalscale << " " << arr[IND(i,j,k)] << std::endl;
    }
    outfile << std::endl;
  }
}


void save_fileXML(const int ns[3], const double ds[3], const double Ls[3], double* arr, char fname[], std::string pathname, double finalscale, double xoffset, double yoffset, std::vector<Electrodes> elec_vec){

  // define major XML nodes
  boost::property_tree::ptree tree;
  boost::property_tree::ptree node_root;       // <sim_out>
  boost::property_tree::ptree node_eng_info;   // <eng_info>
  boost::property_tree::ptree node_sim_params; // <sim_params>
  boost::property_tree::ptree node_electrode;  // <electrode>
  boost::property_tree::ptree node_potential_map;  // <potential>

  // std::string finalpath = fname;
  std::string temp = fname;
  std::string finalpath = pathname+"/"+temp;

  std::cout << finalpath << std::endl;
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
    node_dim.put("<xmlattr>.x1", (std::to_string((elec.x[0]-xoffset)/finalscale/SIMLENGTH).c_str()));
    node_dim.put("<xmlattr>.y1", (std::to_string((elec.y[0]-yoffset)/finalscale/SIMLENGTH).c_str()));
    node_dim.put("<xmlattr>.x2", (std::to_string((elec.x[1]-xoffset)/finalscale/SIMLENGTH).c_str()));
    node_dim.put("<xmlattr>.y2", (std::to_string((elec.y[1]-yoffset)/finalscale/SIMLENGTH).c_str()));
    node_electrode.add_child("dim", node_dim);
    boost::property_tree::ptree node_pot;
    node_pot.put("", std::to_string(elec.potential).c_str());
    node_electrode.add_child("potential", node_pot);
  }

  //potential_map
  const int k = ns[2]/2;
  for (int i = 0; i < ns[0]; i++){
    for (int j = 0; j < ns[1]; j++){
      //create each entry
      boost::property_tree::ptree node_potential_val;
      node_potential_val.put("<xmlattr>.x", (std::to_string((i*ds[0]-xoffset)/finalscale/SIMLENGTH)).c_str());
      node_potential_val.put("<xmlattr>.y", (std::to_string((j*ds[1]-yoffset)/finalscale/SIMLENGTH)).c_str());
      node_potential_val.put("<xmlattr>.val", std::to_string(arr[IND(i,j,k)]).c_str());
      node_potential_map.add_child("potential_val", node_potential_val);
    }
  }
  node_root.add_child("eng_info", node_eng_info);
  // node_root.add_child("sim_params", node_sim_params);
  node_root.add_child("electrode", node_electrode);
  node_root.add_child("potential_map", node_potential_map);
  // node_root.add_child("elec_dist", node_elec_dist);
  tree.add_child("sim_out", node_root);

  // write to file
  boost::property_tree::write_xml(finalpath, tree, std::locale(), boost::property_tree::xml_writer_make_settings<std::string>(' ',4));

  std::cout << "Write to XML complete." << std::endl;

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
          // a[IND(i,j,k)] = PhysConstants::EPS_SI; //Si relative permittivity
          a[IND(i,j,k)] = 1;  //Free space
        } else{
          // a[IND(i,j,k)] = PhysConstants::EPS_SI; //Si relative permittivity
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
        // a[IND(i,j,k)] = 1e16*PhysConstants::QE/PhysConstants::EPS0/eps[IND(i,j,k)]; //in m^-3, scale by permittivity
        a[IND(i,j,k)] = 0; //in m^-3, scale by permittivity
        chi[IND(i,j,k)] = PhysConstants::CHI_SI;
      }
    }
  }
  std::cout << "Finished RHS initialisation" << std::endl;
}
