import numpy as np
import dolfin as df
import matplotlib.pyplot as plt

#General expression
class GenExp(df.Expression):
  def __init__(self, func_in, **kwargs):
    self.data = kwargs["data"]
    self.x = kwargs["x"]
    self.F = func_in
    self.f_vals = []

  def eval(self, val, x):
    f_val = self.F(x) # This should be the thing that goes into the charge density
                      # function. E.g. x[0], or x[1], etc.
    # interpolate the data at the point of evaluation
    val[0] = np.interp(f_val, self.x, self.data)

      # If the expression is vector valued, overload this.
      # def value_shape(self):
      #   return (1,)
