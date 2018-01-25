#ifndef POISSOLVER_H
#define POISSOLVER_H

namespace PhysConstants{
   const double QE = 1.6021766208e-19; //elementary charge
  //  const double PI = 3.14159265358979323846;
   const double EPS0 = 8.85418782e-12; //permittivity of free space
   const double WF_GOLD = 5.1; //workfunction for gold in eV
   const double WF_COPPER = 4.7; //workfunction for copper in eV
   const double WF_ZINC = 4.3; //workfunction for zinc in eV
   const double WF_CESIUM = 2.1; //workfunction for cesium in eV
   const double WF_NICKEL = 5.01; //workfunction for nickel in eV
   const double CHI_SI = 4.05; //electron affinity for silicon in eV from http://www.ioffe.ru/SVA/NSM/Semicond/Si/basic.html
   const double EPS_SI = 11.7; //relative permittivity of silicon
};


class Poissolver{
  // public:
  //   Electrodes(); //default constructor
  //   Electrodes(double, double, double, double, double, double, double, double); //parametrized constructor
  //   ~Electrodes(); //destructor
  //   double x[2]; //xmin (x[0]) and xmax (x[1])
  //   double y[2]; //ymin and ymax
  //   double z[2]; //zmin and zmax
  //   double potential;   //pointer after conversion of vector
  //   double WF;
  //   void draw(const int[3], const double[3], const double[3], double*, std::pair<int,double>*, double*);
};

#endif //POISSOLVER_H
