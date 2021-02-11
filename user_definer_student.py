import numpy as np
# Utility functions to initialize the problem
from Grid.GridProcessing import Grid
from Shapes.ShapesFunctions import *
# Specify the  file that includes dynamic systems
from dynamics.DubinsCapture import *
from dynamics.DubinsCar4D2 import *
# Plot options
from plot_options import *
# Solver core
from solver import HJSolver

import math

""" USER INTERFACES
- Define grid

- Generate initial values for grid using shape functions

- Time length for computations

- Initialize plotting option

- Call HJSolver function
"""


"""
Assign one of the following strings to `compMethod` to specify the characteristics of computation
"none" -> compute Backward Reachable Set
"minVWithV0" -> compute Backward Reachable Tube
"maxVWithVInit" -> compute max V over time
"minVWithVInit" compute min V over time
"""

g = Grid(np.array([-2.0, 0.0, -0.2, -math.pi]), np.array([2.0, 4.0, 1.2, math.pi]), 4, np.array([60, 60, 30, 72]), [3])

# Define my object
my_car = DubinsCar4D2()

# Use the grid to initialize initial value function
Initial_value_f = CylinderShape(g, [3,4], np.array([0, 2, 0, 0]), [1])
Initial_value_f = np.minimum(Initial_value_f, Lower_Half_Space(g, 2, 0))
Initial_value_f = np.minimum(Initial_value_f, Upper_Half_Space(g, 2, 1.1))

# Look-back lenght and time step
lookback_length = 1.0
t_step = 0.05

small_number = 1e-5
tau = np.arange(start=0, stop=lookback_length + small_number, step=t_step)

def solve():
    po = PlotOptions("3d_plot", [0,1,3], [24])
    V = HJSolver(my_car, g, Initial_value_f, tau, "minVWithV0", po)
    # np.save("test_value_fn", V)


if __name__ == "__main__":
   solve()