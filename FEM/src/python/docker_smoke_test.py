"""Check the migrated container through PoisSolver's linear solve methods."""
import math
import matplotlib
matplotlib.use("Agg")
import dolfin as df
from poisson_class import PoissonSolver
from siqadtools import siqadconn

solver = PoissonSolver()
solver.mesh = df.UnitCubeMesh(3, 3, 3)
solver.bc_type = "dirichlet"
solver.eqn = "laplace"
solver.method = "lu"
solver.preconditioner = "default"
solver.init_guess = "zero"
solver.max_abs_error = 1e-12
solver.max_rel_error = 1e-12
solver.max_linear_iters = 1000
solver.setFunctionSpaces()
solver.u_s = df.Function(solver.V)
exact = df.Expression("1+x[0]+2*x[1]+3*x[2]", degree=3)
solver.bcs = [df.DirichletBC(solver.V, exact, "on_boundary")]
solver.F = (
    df.inner(df.grad(solver.u), df.grad(solver.v)) * df.dx
    - df.Constant(0) * solver.v * df.dx
)
solver.setSolverParams()
solver.solve()
error = df.errornorm(exact, solver.u_s, "L2")
if not math.isfinite(error) or error >= 1e-9:
    raise RuntimeError("Linear solve L2 error exceeds tolerance: {}".format(error))
print("DOLFIN {}: PoisSolver linear solve L2 error {}".format(df.__version__, error))
