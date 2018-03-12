// @file:     phys_engine.cc
// @author:   Samuel
// @created:  2017.08.23
// @editted:  2017.08.23 - Samuel
// @license:  GNU LGPL v3
//
// @desc:     Base class for physics engines

#include "phys_connector.h"
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace phys;


//CONSTRUCTOR
PhysicsConnector::PhysicsConnector(const std::string &eng_name_in,
  const std::string &input_path_in, const std::string &output_path_in)
{
  eng_name = eng_name_in;
  input_path = input_path_in;
  output_path = output_path_in;
}


void PhysicsConnector::helloWorld(void)
{
  std::cout << eng_name << ", " << input_path << ", " << output_path << std::endl;
}
