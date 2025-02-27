import sys, os

import math
import numpy as np

from Shapes.RaceTrack import RaceTrack
from Shapes.utils import *
from Grid.GridProcessing import Grid
from dynamics.DubinsCar4D2 import DubinsCar4D2
from plot_options import PlotOptions
from solver import HJSolver

import os, sys, pathlib, json
save_path = os.path.abspath('../safety_sets/')
sys.path.append(save_path)

import pdb

trackName = 'Thruxton'
track = RaceTrack(trackName)
step = 24
interval = 36

# Define my object
my_car = DubinsCar4D2(uMode = 'max')
v_min = 0
v_max = 30

def get_mesh_value(idx, interval=interval):
    start_pt = track.centerline_arr[idx]
    #pdb.set_trace()
    end_idx = (idx + interval) % track.raceline_length 
    end_pt = track.centerline_arr[end_idx]
    yaw = track.race_yaw[idx]
    if idx + interval > track.raceline_length:
        yaw_min = min(np.concatenate([track.race_yaw[idx:], track.race_yaw[:end_idx]])) - math.pi/6
        yaw_max = max(np.concatenate([track.race_yaw[idx:], track.race_yaw[:end_idx]])) + math.pi/6
    else:
        yaw_min = min(track.race_yaw[idx:end_idx]) - math.pi/6
        yaw_max = max(track.race_yaw[idx:end_idx]) + math.pi/6
    
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
            np.array([x_max, y_max+1, v_max, yaw_max]), 4, 
            np.array([15, 15, 15, 18]), [3])

    # The value function should be calculated in the global coordinate.
    basis = np.array([[np.cos(yaw), -np.sin(yaw), start_pt[0]], 
                    [np.sin(yaw), np.cos(yaw), start_pt[1]]])
    V0, XX, YY = track.get_init_value(g, u_init = u_init, basis=basis)
    return g, V0, XX, YY

# Look-back lenght and time step
lookback_length = 0.1 * 6
t_step = 0.1
small_number = 1e-5
tau = np.arange(start=0, stop=lookback_length + small_number, step=t_step)
po = PlotOptions("3d_plot", [0,1,3], [5])

for idx in [7440]:#np.arange(0, track.raceline_length, step = step):
    g, V0, XX, YY = get_mesh_value(idx)
    V1 = HJSolver(my_car, g, V0, tau, "minVWithV0", po, plot_flag = False)
    #pdb.set_trace()
    np.savez_compressed(f"{save_path}/{trackName}/{idx}", 
                        V=V1, 
                        x = g.vs[0][:, 0, 0, 0], 
                        y = g.vs[1][0, :, 0, 0],
                        v = g.vs[2][0, 0, :, 0],
                        yaw = g.vs[3][0, 0, 0, :],
                        xx = XX, yy = YY ## (X, Y) in global coordindate
                        )