#ifndef ELECTRODES_H
#define ELECTRODES_H
#include <utility>

class Electrodes{
  public:
    Electrodes(); //default constructor
    Electrodes(double xmin, double xmax, double ymin, double ymax, double zmin, double zmax, double pot, double wf);
    ~Electrodes(); //destructor
    double x[2]; //xmin (x[0]) and xmax (x[1])
    double y[2]; //ymin and ymax
    double z[2]; //zmin and zmax
    double potential;   //pointer after conversion of vector
    double WF;
    void draw(const int ns[3], const double Ls[3], double* RHS, std::pair<int,double> *electrodemap, double* chi);
    int IND(int i, int j, int k, const int ns[]);
};



#endif //ELECTRODES_H
