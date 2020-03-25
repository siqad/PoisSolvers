import numpy as np
import dolfin as df
import matplotlib.pyplot as plt

class ChargeDensity(df.Expression):
  def __init__(self, func_in, **kwargs):
    self.data = kwargs["data"]
    self.x = kwargs["x"]
    self.F = func_in
    self.f_vals = []
  def eval(self, val, x):
    f_val = self.F(x) # This should be the that goes into the charge density
                      # function. E.g. x[0], or x[1], etc.
    # if f_val not in self.f_vals:
        # self.f_vals.append(f_val)
        # print(self.f_vals)
    # interpolate the data at the point of evaluation
    val[0] = np.interp(f_val, self.x, self.data)

  # If the expression is vector valued, overload this.
  # def value_shape(self):
  #   return (1,)

if __name__ == "__main__":
    rho = np.array([-8.54714149e-07, -8.33732142e-05, -1.60989631e-02, -3.10772450e+00,
     -5.70687012e+02, -1.05369927e+04, -1.10135228e+04, -2.80024598e+02,
      1.07501397e+04,  1.07601003e+04,  8.30053469e+02,  6.00024250e+01,
      4.03255519e+00,  2.23706823e-02,  1.15898807e-04,  6.00232312e-07,
      3.10857530e-09,  1.61008849e-11,  6.83205865e-14,  0.00000000e+00])
    # rho = rho / abs(max(rho))
    # print(rho)
    x = np.array([-1.00000000e-06, -8.94736842e-07, -7.89473684e-07, -6.84210526e-07,
     -5.78947368e-07, -4.73684211e-07, -3.68421053e-07, -2.63157895e-07,
     -1.57894737e-07, -5.26315789e-08,  5.26315789e-08,  1.57894737e-07,
      2.63157895e-07,  3.68421053e-07,  4.73684211e-07,  5.78947368e-07,
      6.84210526e-07,  7.89473684e-07,  8.94736842e-07,  1.00000000e-06])

    # mesh = df.UnitSquareMesh(64, 64)
    print(min(x), max(x))
    mesh = df.RectangleMesh(df.Point(min(x), min(x)), df.Point(max(x), max(x)), 64, 64)
    V = df.FunctionSpace(mesh, 'CG', 1)
    F = df.interpolate(df.Expression('x[0]',  degree=1), V)
    rho_exp = df.interpolate(ChargeDensity(F, data=rho,x=x, degree=1), V)

    df.plot(rho_exp)
    plt.show()
