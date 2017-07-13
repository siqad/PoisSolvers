#ifndef ELECTRODES_H
#define ELECTRODES_H

class Solver;

class Electrodes{
  public:
    Electrodes();
    ~Electrodes();
    double x[2]; //xmin (x[0]) and xmax (x[1])
    double y[2]; //ymin and ymax
    double z[2]; //zmin and zmax
    double potential;   //pointer after conversion of vector
    void draw(const int, const double, const double, double *, std::pair<bool,double> *);

};

#endif //ELECTRODES_H
