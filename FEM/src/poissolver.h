#ifndef POISSOLVER_H
#define POISSOLVER_H

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
        ~PoisSolver(){} //destructor
        void runSolver(void);
        std::string in_path;
        std::string out_path;
  };
}

#endif //POISSOLVER_H
