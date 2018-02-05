// @file:     phys_engine.cc
// @author:   Samuel
// @created:  2017.08.23
// @editted:  2017.08.23 - Samuel
// @license:  GNU LGPL v3
//
// @desc:     Base class for physics engines

#include "phys_engine.h"

#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

using namespace phys;

// CONSTRUCTOR
PhysicsEngine::PhysicsEngine(const std::string &eng_nm, const std::string &i_path, const std::string &o_path)
{
  eng_name = eng_nm;
  problem.readProblem(i_path);
  out_path = o_path;
}

// void PhysicsEngine::writeResultsXml()
// {
//   std::cout << "Write results to XML..." << std::endl;
//   // NOTE in the future, there's probably a range of stuff that can be exported.
//   // for now, only export charge config
//
//   // define major XML nodes
//   bpt::ptree tree;
//   bpt::ptree node_root;       // <sim_out>
//   bpt::ptree node_eng_info;   // <eng_info>
//   bpt::ptree node_sim_params; // <sim_params>
//   bpt::ptree node_physloc;    // <physloc>
//   bpt::ptree node_elec_dist;  // <elec_dist>
//
//   // eng_info
//   node_eng_info.put("engine", eng_name.c_str());
//   node_eng_info.put("version", "TBD"); // TODO real version
//
//   // sim_params
//   // TODO
//
//   // physloc
//   for (auto dbl : db_locs) {
//     bpt::ptree node_dbdot;
//     node_dbdot.put("<xmlattr>.x", std::to_string(dbl.first).c_str());
//     node_dbdot.put("<xmlattr>.y", std::to_string(dbl.second).c_str());
//     node_physloc.add_child("dbdot", node_dbdot);
//   }
//
//   // elec_dist
//   for (auto db_charge : db_charges) {
//     bpt::ptree node_dist;
//     std::string dbc_link;
//     for(auto chg : db_charge)
//       dbc_link.append(std::to_string(chg));
//     node_dist.put("", dbc_link);
//     node_elec_dist.add_child("dist", node_dist);
//   }
//
//   // add nodes to appropriate parent
//   node_root.add_child("eng_info", node_eng_info);
//   node_root.add_child("sim_params", node_sim_params);
//   node_root.add_child("physloc", node_physloc);
//   node_root.add_child("elec_dist", node_elec_dist);
//   tree.add_child("sim_out", node_root);
//
//   // write to file
//   bpt::write_xml(out_path, tree, std::locale(), bpt::xml_writer_make_settings<std::string>(' ',4));
//
//   std::cout << "Write to XML complete." << std::endl;
// }


// void PoisSolver::save_fileXML(double* arr, char fname[], std::vector<Electrodes> elec_vec)
void PhysicsEngine::writeResultsXml()
{
  std::cout << "PhysicsEngine::writeResultsXML()" << std::endl;
  // define major XML nodes
  boost::property_tree::ptree tree;
  boost::property_tree::ptree node_root;       // <sim_out>
  boost::property_tree::ptree node_eng_info;   // <eng_info>
  boost::property_tree::ptree node_sim_params; // <sim_params>
  boost::property_tree::ptree node_electrode;  // <electrode>
  boost::property_tree::ptree node_potential_map;  // <potential>
  // std::string temp = fname;
  // std::string finalpath = SimParams::resultpath + "/" + temp;

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
  boost::property_tree::write_xml(out_path, tree, std::locale(), boost::property_tree::xml_writer_make_settings<std::string>(' ',4));

  std::cout << "Write to XML complete." << std::endl;
}
