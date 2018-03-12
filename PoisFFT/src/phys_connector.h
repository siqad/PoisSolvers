// @file:     phys_engine.h
// @author:   Samuel
// @created:  2017.08.23
// @editted:  2017.08.23 - Samuel
// @license:  GNU LGPL v3
//
// @desc:     Base class for physics engines connectors. Should take care
//            of setting expected problem parameters, parsing problem files,
//            writing result files, etc. Use of the class is recommended, but
//            ultimately optional as users may want to implement their own
//            I/O with the GUI

#ifndef _POIS_SOLVER_PHYS_PHYS_CONNECTOR_H_
#define _POIS_SOLVER_PHYS_PHYS_CONNECTOR_H_

#include "problem.h"
#include "electrodes.h"

#include <string>
#include <vector>
#include <boost/circular_buffer.hpp>

namespace phys{
  namespace bpt = boost::property_tree;
  class PhysicsConnector
  {
  public:
    //CONSTRUCTOR
    PhysicsConnector(const std::string &eng_name_in, const std::string &input_path_in, const std::string &output_path_in);
    //DESTRUCTOR
    ~PhysicsConnector(){};
    void helloWorld(void);
  private:
    std::string eng_name;
    std::string input_path;
    std::string output_path;
  };



}//end namespace phys

#endif
