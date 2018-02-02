#include "poissolver.h"

namespace SimParams
{
  // scaling and offset values
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

int main(int argc,char* argv[]){
  std::cout << "Number of command line arguments: " << argc << std::endl;
  std::vector<Electrodes> elec_vec;
  std::string arg1;
  std::string arg2;
  PoisSolver *ps;

  if(argc == 1){
    std::cout << "No path was passed to the solver program. Program terminating." << std::endl;
  }else if(argc > 2){ //need at LEAST binary call, input path, output path.
    arg1 = argv[1];
    arg2 = argv[2];
    if(arg1.find(".xml") != std::string::npos){ //argv[1] is an xml path, assume it is the INPUT file.
      std::cout << "Input path detected." << std::endl;
      parse_tree(&elec_vec, argv[1]);
      elec_vec = set_buffer(elec_vec);
      std::cout << SimParams::finalscale << " " << SimParams::xoffset << " " << SimParams::yoffset << std::endl;
      if(arg2.find(".xml") != std::string::npos){ //argv[2] is an xml path, assume it is the OUTPUT file.
        std::cout << "Output path detected." << std::endl;
        std::size_t found = arg2.find_last_of("/\\");
        SimParams::resultpath = arg2.substr(0, found);
        std::cout << SimParams::resultpath << std::endl;
        for(int i = 0; i < 1; i++){
          // All the relevant information is inside SimParams.
          ps = new PoisSolver;
          ps->helloWorld();
          ps->worker(i, elec_vec); //where the magic happens
        }
      }else{
        std::cout << "Path not detected. Result XML path needs to be provided as second argument to binary. Terminating." << std::endl;
      }
    }else{
      std::cout << "Path not detected. Problem XML path needs to be provided as first argument to binary. Terminating." << std::endl;
    }
  }
  delete ps;
  return 0;
}

std::vector<Electrodes> set_buffer(std::vector<Electrodes> elec_vec) {
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
  SimParams::finalscale = std::min(xscale, yscale);
  // std::cout << "Final scaling factor is: " << SimParams::finalscale << std::endl;
  for(int i = 0; i < elec_vec.size(); i++){
    elec_vec[i].x[0] *= SimParams::finalscale;
    elec_vec[i].x[1] *= SimParams::finalscale;
    elec_vec[i].y[0] *= SimParams::finalscale;
    elec_vec[i].y[1] *= SimParams::finalscale;
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
  if(xmin < 0.1*SimParams::Ls[0]){  //too far to the left, want positive offset to bring it right
    //find the offset
    SimParams::xoffset = 0.1*SimParams::Ls[0] - xmin;
  }else if(xmax > 0.9*SimParams::Ls[0]){ //too far right, want negative offset to bring it left.
    SimParams::xoffset = 0.9*SimParams::Ls[0] - xmax;
  }
  if(ymin < 0.1*SimParams::Ls[1]){ //too far up
    //find the offset in y
    SimParams::yoffset = 0.1*SimParams::Ls[1] - ymin;
  }else if(ymax > 0.9*SimParams::Ls[1]){ //too far down
    SimParams::yoffset = 0.9*SimParams::Ls[1] - ymax;
  }
  //fix the offsets
  for(int i = 0; i < elec_vec.size(); i++){ //move all points based on offset.
    elec_vec[i].x[0] += SimParams::xoffset;
    elec_vec[i].x[1] += SimParams::xoffset;
    elec_vec[i].y[0] += SimParams::yoffset;
    elec_vec[i].y[1] += SimParams::yoffset;
  }
  return elec_vec;
}

void parse_tree(std::vector<Electrodes> *elecs, std::string path){
  int pix_x1, pix_x2, pix_y1, pix_y2;
  double potential;
  std::cout << "PARSING NOW, PATH NAME IS " << path << std::endl;
  boost::property_tree::ptree tree; // Create empty property tree object
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
