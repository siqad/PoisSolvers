#ifndef ELECTRODES_H
#define ELECTRODES_H
#include <utility>

class Electrodes{
  public:
    Electrodes(); //default constructor
    Electrodes(double, double, double, double, double, double, double, double); //parametrized constructor
    ~Electrodes(); //destructor
    double x[2]; //xmin (x[0]) and xmax (x[1])
    double y[2]; //ymin and ymax
    double z[2]; //zmin and zmax
    double potential;   //pointer after conversion of vector
    double WF;
    void draw(const int[3], const double[3], const double[3], double*, std::pair<int,double>*, double*);
    int IND(int i, int j, int k, const int ns[]);
};



#endif //ELECTRODES_H
