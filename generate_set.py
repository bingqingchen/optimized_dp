import sys, os

import math
import numpy as np

from Shapes.RaceTrack import RaceTrack
from Shapes.utils import *
from Grid.GridProcessing import Grid
from dynamics.DubinsCar4D2 import DubinsCar4D2
from plot_options import PlotOptions
from solver import HJSolver

name = 'Thruxton'
track = RaceTrack(name)
step = 24
interval = 36

# Define my object
my_car = DubinsCar4D2(uMode = 'max')
v_min = 0
v_max = 30

def get_mesh_value(idx, interval=interval):
    start_pt = track.centerline_arr[idx]
    end_pt = track.centerline_arr[idx + interval]
    yaw = track.race_yaw[idx]
    yaw_min = min(track.race_yaw[idx:idx+interval]) - math.pi/6
    yaw_max = max(track.race_yaw[idx:idx+interval]) + math.pi/6
    
    pts = []
    u_init = idx/track.raceline_length
    for pt in [start_pt, end_pt]:
        for tck in [track.tck_in, track.tck_out]:
            _, p_closest, _ = track._calc_shortest_distance(pt, tck, u_init = u_init)
            pts.append(p_closest)

    pts = np.array(pts)
    pts_local = toLocal(pts, start_pt, yaw)
    x_min, x_max, y_min, y_max = fitRectangle(pts_local)

    g = Grid(np.array([x_min, y_min-1, v_min, yaw_min]), 
            np.array([x_max, y_max+1, v_max, yaw_max]), 4, np.array([15, 15, 10, 18]), [3])

    # The value function should be calculated in the global coordinate.
    basis = np.array([[np.cos(yaw), -np.sin(yaw), start_pt[0]], 
                    [np.sin(yaw), np.cos(yaw), start_pt[1]]])
    V0 = track.get_init_value(g, u_init = u_init, basis=basis)
    return g, V0

# Look-back lenght and time step
lookback_length = 0.1 * 6
t_step = 0.1
small_number = 1e-5
tau = np.arange(start=0, stop=lookback_length + small_number, step=t_step)
po = PlotOptions("3d_plot", [0,1,3], [5])

for idx in np.arange(0, track.raceline_length, step = step):
    g, V0 = get_mesh_value(idx)
    HJSolver(my_car, g, V0, tau, "minVWithV0", po)