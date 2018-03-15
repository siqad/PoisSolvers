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

// #include "problem.h"
#include "electrodes.h"

#include <stack>
#include <memory>
#include <map>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
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

    //salutations
    void helloWorld(void);
    //write results to xml, based on the export flags
    void writeResultsXml();

    // Constructor
    void initProblem();

    // File Handling
    bool readProblem();

    // Accessors
    //set a parameter as required for the simulation.
    void setRequiredSimParam(std::string param_name);

    //Input flags
    void setExpectElectrode(bool set_val){expect_electrode = set_val;}
    void setExpectDB(bool set_val){expect_db = set_val;}
    void setExpectAFMPath(bool set_val){expect_afm_path = set_val;}

    //Output flags
    void setExportElecPotential(bool set_val){export_elec_potential = set_val;}
    void setExportDBElecConfig(bool set_val){export_db_elec_config = set_val;}

    //set vector of strings as potential data
    void setElecPotentialData(std::vector<std::vector<std::string>> &data_in){pot_data = data_in;}

    //set vector of strings as electrode data
    void setElectrodeData(std::vector<std::vector<std::string>> &data_in){elec_data = data_in;}

    //get the required simulation parameter vector.
    std::vector<std::string> getRequiredSimParam(void){return req_params;}

    //! Checks if a parameter exists given the parameter key.
    bool parameterExists(const std::string &key) {return sim_params.find(key) != sim_params.end();}

    //! Getter for a parameter, given a parameter key.
    std::string getParameter(const std::string &key) {return sim_params.find(key) != sim_params.end() ? sim_params.at(key) : "";}

    //simulation inputs and outputs
    std::vector<std::vector<std::string>> pot_data;
    std::vector<std::vector<std::string>> elec_data;

    std::vector<std::pair<float,float>> db_locs;
    boost::circular_buffer<std::vector<int>> db_charges;

    // STRUCTS
    // electrode
    struct Electrode {
      float x1,x2,y1,y2;      // pixel location of electrode.
      float potential;  // voltage that the electrode is set to
      Electrode(float in_x1, float in_x2, float in_y1, float in_y2, float in_potential)
        : x1(in_x1), x2(in_x2), y1(in_y1), y2(in_y2), potential(in_potential) {};
    };

    // aggregate
    class Aggregate
    {
    public:
      std::vector<std::shared_ptr<Aggregate>> aggs;
      std::vector<std::shared_ptr<Electrode>> elecs;

      // Properties
      int size(); // returns the number of contained elecs, including those in children aggs
    };

    // ITERATOR
    typedef std::vector<std::shared_ptr<Electrode>>::const_iterator ElecIter;
    typedef std::vector<std::shared_ptr<Aggregate>>::const_iterator AggIter;

    class ElecIterator
    {
    public:
      explicit ElecIterator(std::shared_ptr<Aggregate> root, bool begin=true);

      //~Iterator() {delete agg_stack;};

      ElecIterator& operator++(); // recursive part here
      bool operator==(const ElecIterator &other) {return other.elec_iter == elec_iter;}
      bool operator!=(const ElecIterator &other) {return other.elec_iter != elec_iter;}
      std::shared_ptr<Electrode> operator*() const {return *elec_iter;}
    private:
      ElecIter elec_iter;                   // points to the current electrode
      std::shared_ptr<Aggregate> curr;  // current working Aggregate
      std::stack<std::pair<std::shared_ptr<Aggregate>, AggIter>> agg_stack;
      // add a new aggregate pair to the stack
      void push(std::shared_ptr<Aggregate> agg);
      // pop the aggregate stack
      void pop();
    };
    ElecIterator begin() {return ElecIterator(elec_tree);}
    ElecIterator end() {return ElecIterator(elec_tree, false);}

  private:

    bool readProgramProp(const bpt::ptree &);
    bool readMaterialProp(const bpt::ptree &);
    bool readSimulationParam(const bpt::ptree &sim_params_tree);
    bool readDesign(const bpt::ptree &subtree, const std::shared_ptr<Aggregate> &agg_parent);
    bool readItemTree(const bpt::ptree &subtree, const std::shared_ptr<Aggregate> &agg_parent);
    bool readElectrode(const bpt::ptree &subtree, const std::shared_ptr<Aggregate> &agg_parent);

    // Variables
    std::shared_ptr<Aggregate> elec_tree;
    std::map<std::string, std::string> program_props;
    // std::map<std::string, std::string> material_props; TODO probably need a different structure for this
    std::map<std::string, std::string> sim_params;


    ElecIter elec_iter;                   // points to the current electrode
    std::shared_ptr<Aggregate> curr;  // current working Aggregate
    std::stack<std::pair<std::shared_ptr<Aggregate>, AggIter>> agg_stack;


    std::string eng_name;
    std::string input_path;
    std::string output_path;
    std::vector<std::string> req_params;

    bool expect_electrode;
    bool expect_db;
    bool expect_afm_path;

    bool export_elec_potential;
    bool export_db_elec_config;
  };



}//end namespace phys

#endif
