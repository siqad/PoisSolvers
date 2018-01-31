#include "poissolver.h"

int main(int argc,char* argv[]){
  std::cout << "Number of command line arguments: " << argc << std::endl;
  std::vector<Electrodes> elec_vec;
  std::string arg1;
  std::string arg2;
  std::string outpath;
  double finalscale = 0;
  double xoffset = 0;
  double yoffset = 0;

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
        for(int i = 0; i < 1; i++){
          worker(i, elec_vec, outpath, finalscale, xoffset, yoffset); //where the magic happens
        }
      }else{
        std::cout << "Path not detected. Result XML path needs to be provided as second argument to binary. Terminating." << std::endl;
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
  if(xmin < 0.1*SimParams::Ls[0] || xmax > 0.9*SimParams::Ls[0]){
    xlength = xmax - xmin;
    xscale = 0.8*SimParams::Ls[0]/xlength;
    std::cout << "xlength: " << xlength << std::endl;
  }
  if(ymin < 0.1*SimParams::Ls[1] || ymax > 0.9*SimParams::Ls[1]){
    ylength = ymax - ymin;
    yscale = 0.8*SimParams::Ls[1]/ylength;
    std::cout << "ylength: " << ylength << std::endl;
  }
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

  if(xmin < 0.1*SimParams::Ls[0]){  //too far to the left, want positive offset to bring it right
    //find the offset
    *xoffset = 0.1*SimParams::Ls[0] - xmin;
  }else if(xmax > 0.9*SimParams::Ls[0]){ //too far right, want negative offset to bring it left.
    *xoffset = 0.9*SimParams::Ls[0] - xmax;
  }
  if(ymin < 0.1*SimParams::Ls[1]){ //too far up
    //find the offset in y
    *yoffset = 0.1*SimParams::Ls[1] - ymin;
  }else if(ymax > 0.9*SimParams::Ls[1]){ //too far down
    *yoffset = 0.9*SimParams::Ls[1] - ymax;
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
                elecs->push_back(Electrodes(pix_x1*SimParams::Ls[0], pix_x2*SimParams::Ls[0], pix_y1*SimParams::Ls[1], pix_y2*SimParams::Ls[1], 0.3e-6, 0.7e-6, potential, PhysConstants::WF_GOLD));
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

  save_file2D(arr, SimParams::OUTFILE, resultpath, finalscale, xoffset, yoffset); //solution is in arr
  save_file2D(RHS, SimParams::RHOFILE, resultpath, finalscale, xoffset, yoffset);
  save_file2D(correction, SimParams::CORRFILE, resultpath, finalscale, xoffset, yoffset);
  save_fileXML(arr, SimParams::RESXML, resultpath, finalscale, xoffset, yoffset, elec_vec);

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

void calc_charge(double* RHS , std::vector<Electrodes> elecs){
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

void create_electrode(double* RHS, std::pair<int,double> *electrodemap, double* chi, std::vector<Electrodes> elecs){
  for(int i = 0; i < elecs.size(); i++){
    elecs[i].draw(SimParams::ns, SimParams::ds, SimParams::Ls, RHS, electrodemap, chi); //separately call draw for each electrode.
  }
}

void apply_correction(double *RHS, double *correction, std::pair<int,double> *electrodemap){
  for(int i = 0; i < SimParams::ns[0]*SimParams::ns[1]*SimParams::ns[2]; i++){
    if(electrodemap[i].first == true){ //only correct error at electrode surfaces.
      RHS[i] -= correction[i];
    }
  }
}

double check_error(double *arr, double *correction, std::pair<int,double> *electrodemap, int *indexErr, double *eps){
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

void save_file2D(double* arr, char fname[], std::string pathname, double finalscale, double xoffset, double yoffset){
  // std::string finalpath = fname;
  std::string temp = fname;
  std::string finalpath = pathname+"/"+temp;
  std::ofstream outfile;
  // std::cout << pathname << std::endl;
  // std::cout << finalpath << std::endl;
  // std::cout << temp << std::endl;
  outfile.open(finalpath.c_str(), std::ios_base::out | std::ios_base::trunc );
  std::cout << "Dumping to " << finalpath << std::endl;
  const int k = SimParams::ns[2]/2;
  for (int i = 0; i < SimParams::ns[0]; i++){
    for (int j = 0; j < SimParams::ns[1]; j++){
      outfile << std::setprecision(5) << std::scientific << (i*SimParams::ds[0]-xoffset)/finalscale << " "
              << (j*SimParams::ds[1]-yoffset)/finalscale << " " << arr[SimParams::IND(i,j,k)] << std::endl;
    }
    outfile << std::endl;
  }
}


void save_fileXML(double* arr, char fname[], std::string pathname, double finalscale, double xoffset, double yoffset, std::vector<Electrodes> elec_vec){

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
    node_dim.put("<xmlattr>.x1", (std::to_string((elec.x[0]-xoffset)/finalscale/SimParams::Ls[0]).c_str()));
    node_dim.put("<xmlattr>.y1", (std::to_string((elec.y[0]-yoffset)/finalscale/SimParams::Ls[1]).c_str()));
    node_dim.put("<xmlattr>.x2", (std::to_string((elec.x[1]-xoffset)/finalscale/SimParams::Ls[0]).c_str()));
    node_dim.put("<xmlattr>.y2", (std::to_string((elec.y[1]-yoffset)/finalscale/SimParams::Ls[1]).c_str()));
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
      node_potential_val.put("<xmlattr>.x", (std::to_string((i*SimParams::ds[0]-xoffset)/finalscale/SimParams::Ls[0])).c_str());
      node_potential_val.put("<xmlattr>.y", (std::to_string((j*SimParams::ds[1]-yoffset)/finalscale/SimParams::Ls[1])).c_str());
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


void init_eps(double* eps){
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
  std::cout << "Finished rho initialisation" << std::endl;
}

void init_rhs(double* chi, double* eps, double* rhs){
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
