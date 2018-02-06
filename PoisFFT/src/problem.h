// @file:     problem.h
// @author:   Samuel
// @created:  2017.08.23
// @editted:  2017.08.23 - Samuel
// @license:  GNU LGPL v3
//
// @desc:     Definition of the problem - dbdot loc, material properties, etc.

#ifndef _SIM_ANNEAL_PHYS_PROBLEM_H_
#define _SIM_ANNEAL_PHYS_PROBLEM_H_

#include <vector>
#include <stack>
#include <memory>
#include <string>
#include <map>
#include <iostream>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace phys{

  namespace bpt = boost::property_tree;

  class Problem
  {
  public:
    // Constructor
    Problem(const std::string &fname);
    Problem() {initProblem();};
    void initProblem();

    // Destructor
    ~Problem() {};

    // File Handling
    bool readProblem(const std::string &fname);

    // Accessors
    bool parameterExists(const std::string &key) {return sim_params.find(key) != sim_params.end();}
    std::string getParameter(const std::string &key) {return sim_params.find(key) != sim_params.end() ? sim_params.at(key) : "";}

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
  };
}

#endif
