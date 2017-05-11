#include <dolfin.h>
#include "Poisson.h"
//#include <cmath>

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
    public:
    double zz, s0, sfluc, eeps;
    Source(double z, double si, double sf):Expression() {
        zz=z; s0=si; sfluc=sf; eeps=1e-20;
    }

  void eval(Array<double>& values, const Array<double>& x) const
  {
      if (abs(x[2]-zz)<1e-9) {
          double rv=-1.0+2.0*dolfin::rand();
          values[0]=s0+rv*sfluc;
      }
      else {
          values[0]=0;
      }
  }
};

// Boundary flux (Neumann boundary condition)
class Flux : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    values[0] = 0;
  }
};
// (Dirichlet boundary condition)
class DirichletBoundary : public SubDomain {
    public:
        double zz, vol;
        DirichletBoundary(double z, double v) : SubDomain() {
            zz=z; vol=v;
        }
    bool inside(const Array<double>& x) const {
//  bool inside(const Array<double>& x, bool on_boundary) const {
        return (x[2]-zz>=0) && (x[2]-zz<2e-9);
    }
};

int main()
{
    double epsilon=8.854e-12;
    double Lx=500e-9;
    double Ly=500e-9;
    double Lz=20e-9;
    double zTI=10e-9;
    double zele=10e-9;
    double Nx=20;
    double Ny=20;
    double Nz=20;
    double s0=0;
    double sfluc=1e13;
    Mesh mesh;
  // Create mesh and function space
  //BoxMesh mesh(0,0,0,Lx,Ly,Lz,Nx,Ny,Nz);
    mesh = BoxMesh(Point(0,0,0),Point(10,10,10), 10, 20, 30);
  //BoxMesh mesh(0,0,0,20,20,20,20,20,20);
  //plot(mesh);
  //interactive();
  //std::getchar();
  Poisson::FunctionSpace V(mesh);
  //V.print_dofmap();
  cout << V.dim() << endl;
  /* BC */
  Constant v0(0.0,0.0);
  DirichletBoundary dBC(zele,1);
  DirichletBC dirbc(V,v0,dBC);

  // Define variational problem
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);
  Source f(zTI,s0,sfluc/epsilon);
  Flux g;
  L.f = f;
  L.g = g;
  Function w(V);
  //LinearVariationalProblem problem(a, L, w);

  // Compute solution
  //LinearVariationalSolver solver(problem);
  //solver.solve();

    solve(a==L, w, dirbc);

  // Extract subfunction
  Function u = w[0];

  // Plot solution
  plot(u);
  interactive();

  return 0;
}
