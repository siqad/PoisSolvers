#ifndef POISSOLVER_H
#define POISSOLVER_H



#include "phys_connector.h"
#include "electrodes.h" // Electrodes
#include <cmath> // fabs()
#include <iostream> // cout
#include <utility> // pair
#include <algorithm> // min(), max()
#include <ctime> // clock()
#include <boost/property_tree/ptree.hpp> // XML stuff
#include <boost/property_tree/xml_parser.hpp> // XML stuff
#include <vector> // vector
#include <string> // string
#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>

namespace phys{
  class PoisSolver
  {
    public:
        PoisSolver(const std::string& i_path, const std::string& o_path);
        ~PoisSolver(){} //destructor1
        void initSolver(void);
        void runSolver(void);
        void loadPotentialTxt(std::vector<std::vector<float>> &elec_data);
        void initVars(void);
        void exportData(void);        
        std::vector<Electrodes> elec_vec; // location of elecs
        std::string bc;
        int resolution;
        double length;
        double max_error;
        double z_offset;
        double z_thickness;
        std::string mode;
        double high_pot;
        double low_pot;
        int steps;
        std::vector<std::vector<float>> e_pot_2D;
        PhysicsConnector* phys_con;
  };
}

double interpolate(double q11, double q12, double q21,  double q22, double x1, double x2, double y1, double y2, double x, double y);
#endif //POISSOLVER_H
