/*
    poisson.cpp (poisson)
    Jess Robertson, 2010-11-17

    Write methods for Poisson class
*/

#include "poisson.hpp"

using namespace mgrid;

Poisson::Poisson(const Settings& settings): LinearMultigrid::LinearMultigrid(settings) {
    nxfine = solution[finestLevel].rows();
    nzfine = solution[finestLevel].columns();
    std::cout << " -- Initial dimensions: (" << nxfine << ", " << nzfine
              << "), with aspect: " << settings.aspectRatio << std::endl;

    // Set boundary conditions for velocity array
    //    solution.boundaryConditions.set(leftBoundary,   zeroNeumannCondition);
    solution.boundaryConditions.set(leftBoundary,   zeroDirichletCondition);
    solution.boundaryConditions.set(rightBoundary,  zeroDirichletCondition);
    //    solution.boundaryConditions.set(topBoundary,    zeroNeumannCondition);
    solution.boundaryConditions.set(topBoundary,    zeroDirichletCondition);
    solution.boundaryConditions.set(bottomBoundary, zeroDirichletCondition);

    // Set source
      source[finestLevel] = 12.00; sourceIsSet = true;
      std::cout << source[finestLevel] << std::endl;
      std::cout << "Right here" << std::endl;
}
Poisson::~Poisson() { /* pass */ }

std::string Poisson::filename(std::string root) {

    std::ofstream outfile;
    outfile.open("test.txt", std::ios_base::out | std::ios_base::trunc );
    std::ostringstream name;
    name.precision(1);  // Print variables to one decimal place
    name << root << "A" << std::fixed << aspect;
    //std::cout << std::setprecision(5) << std::fixed << solution[finestLevel] << std::endl;
    outfile << std::setprecision(5) << std::fixed << solution[finestLevel] << std::endl;
    return name.str();
}

void Poisson::solve() {
    // Solve using linear multigrid method
    multigrid();
}
