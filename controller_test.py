from user_definer_answer import g as g_answer
from user_definer_answer import my_car as car_answer
from user_definer_answer import Initial_value_f as ivf_answer
from user_definer_answer import tau as tau_answer

from user_definer_student import g as g_student
from user_definer_student import my_car as car_student
from user_definer_student import Initial_value_f as ivf_student
from user_definer_student import tau as tau_student


from solver import HJSolver
from plot_options import PlotOptions
from Plots.plotting_utilities import plot_isosurface
import numpy as np
import sys
from os import path


# The following spa_derivX*_4d are from computeGraphs/graph_4D.py
# TODO: figure how to use functions with heterocl code without having to build 
# #     them, e.g. hcl.if_()...
def spa_derivX1_4d(i, j, k, l, V, g):  # Left -> right == Outer Most -> Inner Most
    left_deriv = 0.0
    right_deriv = 0.0
    if 0 not in g.pDim:
        if i == 0:
            left_boundary = 0.0
            left_boundary = V[i, j, k, l] + abs(
                V[i + 1, j, k, l] - V[i, j, k, l]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[0]
            right_deriv = (V[i + 1, j, k, l] - V[i, j, k, l]) / g.dx[0]
        elif i == V.shape[0] - 1:
            right_boundary = 0.0
            right_boundary = V[i, j, k, l] + abs(
                V[i, j, k, l] - V[i - 1, j, k, l]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - V[i - 1, j, k, l]) / g.dx[0]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[0]
        elif i != 0 and i != V.shape[0] - 1:
            left_deriv = (V[i, j, k, l] - V[i - 1, j, k, l]) / g.dx[0]
            right_deriv = (V[i + 1, j, k, l] - V[i, j, k, l]) / g.dx[0]
        return left_deriv, right_deriv
    else:
        if i == 0:
            left_boundary = 0.0
            left_boundary = V[V.shape[0] - 1, j, k, l]
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[0]
            right_deriv = (V[i + 1, j, k, l] - V[i, j, k, l]) / g.dx[0]
        elif i == V.shape[0] - 1:
            right_boundary = 0.0
            right_boundary = V[0, j, k, l]
            left_deriv = (V[i, j, k, l] - V[i - 1, j, k, l]) / g.dx[0]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[0]
        elif i != 0 and i != V.shape[0] - 1:
            left_deriv = (V[i, j, k, l] - V[i - 1, j, k, l]) / g.dx[0]
            right_deriv = (V[i + 1, j, k, l] - V[i, j, k, l]) / g.dx[0]
        return left_deriv, right_deriv

def spa_derivX2_4d(i, j, k, l, V, g): # Left -> right == Outer Most -> Inner Most
    left_deriv = 0.0
    right_deriv = 0.0
    if 1 not in g.pDim:
        if j == 0:
            left_boundary = 0.0
            left_boundary = V[i, j, k, l] + abs(
                V[i, j + 1, k, l] - V[i, j, k, l]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[1]
            right_deriv = (V[i, j + 1, k, l] - V[i, j, k, l]) / g.dx[1]
        elif j == V.shape[1] - 1:
            right_boundary = 0.0
            right_boundary = V[i, j, k, l] + abs(
                V[i, j, k, l] - V[i, j - 1, k, l]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - V[i, j - 1, k, l]) / g.dx[1]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[1]
        elif j != 0 and j != V.shape[1] - 1:
            left_deriv = (V[i, j, k, l] - V[i, j - 1, k, l]) / g.dx[1]
            right_deriv = (V[i, j + 1, k, l] - V[i, j, k, l]) / g.dx[1]
        return left_deriv, right_deriv
    else:
        if j == 0:
            left_boundary = 0.0
            left_boundary = V[i, V.shape[1] - 1, k, l]
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[1]
            right_deriv = (V[i, j + 1, k, l] - V[i, j, k, l]) / g.dx[1]
        elif j == V.shape[1] - 1:
            right_boundary = 0.0
            right_boundary = V[i, 0, k, l]
            left_deriv = (V[i, j, k, l] - V[i, j - 1, k, l]) / g.dx[1]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[1]
        elif j != 0 and j != V.shape[1] - 1:
            left_deriv = (V[i, j, k, l] - V[i, j - 1, k, l]) / g.dx[1]
            right_deriv = (V[i, j + 1, k, l] - V[i, j, k, l]) / g.dx[1]
        return left_deriv, right_deriv




def spa_derivX3_4d(i, j, k, l, V, g):  # Left -> right == Outer Most -> Inner Most
    left_deriv = 0.0
    right_deriv = 0.0
    if 2 not in g.pDim:
        if k == 0:
            left_boundary = 0.0
            left_boundary = V[i, j, k, l] + abs(
                V[i, j, k + 1, l] - V[i, j, k, l]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[2]
            right_deriv = (V[i, j, k + 1, l] - V[i, j, k, l]) / g.dx[2]
        elif k == V.shape[2] - 1:
            right_boundary = 0.0
            right_boundary = V[i, j, k, l] + abs(
                V[i, j, k, l] - V[i, j, k - 1, l]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - V[i, j, k - 1, l]) / g.dx[2]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[2]
        elif k != 0 and k != V.shape[2] - 1:
            left_deriv = (V[i, j, k, l] - V[i, j, k - 1, l]) / g.dx[2]
            right_deriv = (V[i, j, k + 1, l] - V[i, j, k, l]) / g.dx[2]
        return left_deriv, right_deriv
    else:
        if k == 0:
            left_boundary = 0.0
            left_boundary = V[i, j, V.shape[2] - 1, l]
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[2]
            right_deriv = (V[i, j, k + 1, l] - V[i, j, k, l]) / g.dx[2]
        elif k == V.shape[2] - 1:
            right_boundary = 0.0
            right_boundary = V[i, j, 0, l]
            left_deriv = (V[i, j, k, l] - V[i, j, k - 1, l]) / g.dx[2]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[2]
        elif k != 0 and k != V.shape[2] - 1:
            left_deriv = (V[i, j, k, l] - V[i, j, k - 1, l]) / g.dx[2]
            right_deriv = (V[i, j, k + 1, l] - V[i, j, k, l]) / g.dx[2]
        return left_deriv, right_deriv

def spa_derivX4_4d(i, j, k, l, V, g):  # Left -> right == Outer Most -> Inner Most
    left_deriv = 0.0
    right_deriv = 0.0
    if 3 not in g.pDim:
        if l == 0:
            left_boundary = 0.0
            left_boundary = V[i, j, k, l] + abs(
                V[i, j, k, l + 1] - V[i, j, k, l]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[3]
            right_deriv = (V[i, j, k, l + 1] - V[i, j, k, l]) / g.dx[3]
        elif l == V.shape[3] - 1:
            right_boundary = 0.0
            right_boundary = V[i, j, k, l] + abs(
                V[i, j, k, l] - V[i, j, k, l - 1]
            ) * np.sign(V[i, j, k, l])
            left_deriv = (V[i, j, k, l] - V[i, j, k, l - 1]) / g.dx[3]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[3]
        elif l != 0 and l != V.shape[3] - 1:
            left_deriv = (V[i, j, k, l] - V[i, j, k, l - 1]) / g.dx[3]
            right_deriv = (V[i, j, k, l + 1] - V[i, j, k, l]) / g.dx[3]
        return left_deriv, right_deriv
    else:
        if l == 0:
            left_boundary = 0.0
            left_boundary = V[i, j, k, V.shape[3] - 1]
            left_deriv = (V[i, j, k, l] - left_boundary) / g.dx[3]
            right_deriv = (V[i, j, k, l + 1] - V[i, j, k, l]) / g.dx[3]
        elif l == V.shape[3] - 1:
            right_boundary = 0.0
            right_boundary = V[i, j, k, 0]
            left_deriv = (V[i, j, k, l] - V[i, j, k, l - 1]) / g.dx[3]
            right_deriv = (right_boundary - V[i, j, k, l]) / g.dx[3]
        elif l != 0 and l != V.shape[3] - 1:
            left_deriv = (V[i, j, k, l] - V[i, j, k, l - 1]) / g.dx[3]
            right_deriv = (V[i, j, k, l + 1] - V[i, j, k, l]) / g.dx[3]
        return left_deriv, right_deriv


def calculate_optimal_control(state, V, grid, car):
    """
    recall: to calculate optimal control we need the following
        u^* = argmax_{u} (∂V/∂x)^T • f(x, u)
            = argmax_{u} [∂V/∂x_1] • f(x, u) 
                         [∂V/∂x_2]
                         [∂V/∂x_3]
                         [∂V/∂x_4]
        where V := value function from HJSolver
              x := state
              u := control input
              f := dynamic system
    """

    i, j, k, l = state

    # calculate (∂V/∂x)^T
    dV_dx1_L, dV_dx1_R = spa_derivX3_4d(i, j, k, l, V, grid)
    dV_dx2_L, dV_dx2_R = spa_derivX3_4d(i, j, k, l, V, grid)
    dV_dx3_L, dV_dx3_R = spa_derivX3_4d(i, j, k, l, V, grid)
    dV_dx4_L, dV_dx4_R = spa_derivX4_4d(i, j, k, l, V, grid)

    dV_dx1 = (dV_dx1_L + dV_dx1_R) / 2
    dV_dx2 = (dV_dx2_L + dV_dx2_R) / 2
    dV_dx3 = (dV_dx3_L + dV_dx3_R) / 2
    dV_dx4 = (dV_dx4_L + dV_dx4_R) / 2


    # calculate u*
    _, _, opt_a, opt_w = car.opt_ctrl(0, state, (dV_dx1, dV_dx2, dV_dx3, dV_dx4), use_hcl = False)

    return opt_a, opt_w


def compare_V(V_correct, V_test):
    """ Compares two Value functions of same size

    Args:
        V_correct (np.array): The correct value function
        V_test (np.array): The testing value function
    """
    # the following plots show the minimal backwards reacable tube
    # the green shape represents the boundary of it
    # anything within the green shapes represents states that are
    # guarentied to collide with an object no matter how well we control the car
    # anything outside of the green shape represents states that will never collide
    # with an object.

    # since our grid is 4D, we fix the second dimention to constant speed
    # try changing the value to see how the BRT changes as the car increase/decreases in 
    # speed
    # po = PlotOptions("3d_plot", [0,1,3], [24])
    # plot_isosurface(g, V, po)

    # the following plot shows the constrained speed from the half spaces
    # po = PlotOptions("3d_plot", [0,1,2], [40])
    # plot_isosurface(g, V, po)


    # assume that the grid was the same for both Value functions
    total_pts = np.prod(g_answer.pts_each_dim)
    epsilon = 0.001
    correct_opt_ctrl = 0

    # for i in range(g_answer.pts_each_dim[0]):
    #     for j in range(g_answer.pts_each_dim[1]):
    #         for k in range(g_answer.pts_each_dim[2]):
    #             for l in range(g_answer.pts_each_dim[3]):
    #                 state = (i, j, k, l)
    #                 opt_a_correct, opt_w_correct = calculate_optimal_control(state, V_correct, g_answer, car_answer) 
    #                 opt_a, opt_w = calculate_optimal_control(state, V_test, g_answer, car_answer)
    #                 if opt_a_correct == opt_a and opt_w_correct == opt_w:
    #                     correct_opt_ctrl += 1
    # print(f"# of correct optimal control: {(correct_opt_ctrl / total_pts)*100:2f}")

    correct_value = (abs(V_correct - V_test) < epsilon).sum()
    print(f"# of correct value : {(correct_value / total_pts)*100:2f}")




if __name__ == "__main__":
    # correct value function
    po = PlotOptions("3d_plot", [0,1,3], [24])
    print("calculating true value function")
    V1 = HJSolver(car_answer, g_answer, ivf_answer, tau_answer, "minVWithV0", po)

    print("calculating other value function")
    V2 = HJSolver(car_student, g_student, ivf_student, tau_student, "minVWithV0", po)
    # assert(V1.shape == V2.shape)
    compare_V(V1, V2)