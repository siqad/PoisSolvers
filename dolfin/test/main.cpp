// Copyright (C) 2006-2011 Anders Logg
//
// This file is part of DOLFIN.
//
// DOLFIN is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// DOLFIN is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2006-02-07
// Last changed: 2013-03-11
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = 10*exp(-((x - 0.5)^2 + (y - 0.5)^2) / 0.02)
//
// and boundary conditions given by
//
//     u(x, y) = 0        for x = 0 or x = 1
// du/dn(x, y) = sin(5*x) for y = 0 or y = 1

#include <dolfin.h>
#include "Poisson.h"

using namespace dolfin;

// Source term (right-hand side)
class Source : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    double dx = x[0] - 0.5;
    double dy = x[1] - 0.5;
    double dz = x[2] - 0.5;
    values[0] = 10*exp(-(dx*dx + dy*dy + dz*dz) / 0.02);
  }
};

// Normal derivative (Neumann boundary condition)
class dUdN : public Expression
{
  void eval(Array<double>& values, const Array<double>& x) const
  {
    //values[0] = sin(5*x[0]);
    values[0] = 100*x[0];
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const Array<double>& x, bool on_boundary) const
  {
    return x[0] < DOLFIN_EPS or x[0] > 1.0 - DOLFIN_EPS;
  }
};

int main()
{
  // Create mesh and function space
  //UnitSquareMesh mesh(32, 32);
  int nx = 10;
  int ny = 20;
  int nz = 30;
  Point p0(0, 0, 0);
  Point p1(1, 2, 3);
  BoxMesh mesh(p0,p1, nx, ny, nz);
  Poisson::FunctionSpace V(mesh);

  // Define boundary condition
  Constant u0(0.0, 0.0);
  Constant u1(10, 20);
  DirichletBoundary boundary0;
  DirichletBoundary boundary1;
  DirichletBC bc0(V, u0, boundary0);
  DirichletBC bc1(V, u1, boundary1);
  // Define variational forms
  Poisson::BilinearForm a(V, V);
  Poisson::LinearForm L(V);

  Source f;
  dUdN g;
  L.f = f;
  L.g = g;

  // Compute solution
  Function u(V);
  solve(a == L, u, bc0, bc1);

  //extract potential subfunction
  Function w = u[0];
  // Save solution in VTK format
  File file("poisson.pvd");
  //file << u;
  file << w;

  // Plot solution
  //plot(u);
  plot(w);
  interactive();

  return 0;
}
