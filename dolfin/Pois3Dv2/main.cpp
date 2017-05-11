#include <dolfin.h>
//#include "RT0.h"

using namespace dolfin;

// Boundaries
class SE : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[1],0.0,DOLFIN_EPS) && x[0] >= 1;}
};

class SW : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[0],0.0,DOLFIN_EPS) && x[1] >= 1;}
};

class NW : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[1],100.0,DOLFIN_EPS);}
};

class NE : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[2],100.0,DOLFIN_EPS);}
};

class Top : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[2], 30.0, DOLFIN_EPS) && (x[0] <= 99 || x[1] <= 99);}
};

class Bottom: public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[2], 0.0, DOLFIN_EPS);}
};

class SEpressure : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[1],0.0,DOLFIN_EPS) && x[0] <= 1;}
};

class SWpressure : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[0],0.0,DOLFIN_EPS) && x[1] <= 1;}
};

class Toppressure : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  { return near(x[2], 30.0, DOLFIN_EPS) && (x[0] >= 99 && x[1] >= 99);}
};

// Boundary source for flux boundary condition
class BoundarySource : public Expression
{
  public:

    BoundarySource(const Mesh& mesh) : Expression(3), mesh(mesh) {}

    void eval(Array<double>& values, const Array<double>& x,
              const ufc::cell& ufc_cell) const
    {
      dolfin_assert(ufc_cell.local_facet >= 0);

      Cell cell(mesh, ufc_cell.index);
      Point n = cell.normal(ufc_cell.local_facet);

      const double g = 0;
      values[0] = g*n[0];
      values[1] = g*n[1];
      values[2] = g*n[2];
    }

  private:
    const Mesh& mesh;
};

int main(int argc, char* argv[])
{
  // Read any command line flags
  parameters.parse(argc,argv);

  // Create mesh
  BoxMesh mesh(0,0,0,100,100,30,200,200,5);

  // Construct function space
  W = FunctionSpace::FunctionSpace(mesh);
  a = BilinearForm(W, W);
  L = LinearForm(W);

  // Define boundaries
  SubSpace W0(W, 0);
  BoundarySource G(mesh);
  SW sw;
  SE se;
  NE ne;
  NW nw;
  Top top;
  Bottom bottom;
  SWpressure swpressure;
  SEpressure sepressure;
  Toppressure toppressure;
  MeshFunction<std::size_t> boundaries(mesh,mesh.topology().dim()-1);
  boundaries.set_all(0);
  sw.mark(boundaries,1);
  se.mark(boundaries,2);
  ne.mark(boundaries,3);
  nw.mark(boundaries,4);
  top.mark(boundaries,5);
  bottom.mark(boundaries,6);
  swpressure.mark(boundaries,7);
  sepressure.mark(boundaries,8);
  toppressure.mark(boundaries,9);

  // Apply boundary conditions
  DirichletBC bc1(W0,G,sw);
  DirichletBC bc2(W0,G,se);
  DirichletBC bc3(W0,G,ne);
  DirichletBC bc4(W0,G,nw);
  DirichletBC bc5(W0,G,top);
  DirichletBC bc6(W0,G,bottom);
  std::vector<const DirichletBC*> bcs;
  bcs.push_back(&bc1);
  bcs.push_back(&bc2);
  bcs.push_back(&bc3);
  bcs.push_back(&bc4);
  bcs.push_back(&bc5);
  bcs.push_back(&bc6);

  // Material properties
  Constant alpha(3.95e7);
  Constant rhob(0.0,0.0,-4699);
  Constant p_Production(101325.0);
  Constant p_Injection(101325000.0);
  a.alpha = alpha;
  L.rhob = rhob;
  L.p_Production = p_Production;
  L.p_Injection = p_Injection;

  // Compute solution
  Function w(W);
  Matrix A;
  Vector b;
  assemble_system(A,b,a,L,bcs);
  solve(A,*w.vector(),b,"gmres","bjacobi");

  // Extract sub functions (function views)
  Function& sigma = w[0];
  Function& u = w[1];
  list_timings();
  // Plot solutions
  //plot(u);
  //plot(sigma);
  //interactive();

  return 0;
}
