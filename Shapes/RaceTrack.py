import math
import numpy as np
import scipy.interpolate as si
import scipy.optimize as so, scipy.spatial.distance as ssd, scipy.integrate

import os, sys, pathlib, json
l2r_path = os.path.abspath('../learn-to-race')
sys.path.append(l2r_path)

from Shapes.utils import *

class RaceTrack():
    def __init__(self, trackName):
        # Map the track name to the location of the file
        if trackName.lower() == 'thruxton':
            trackFile = os.path.join(l2r_path, 'l2r/racetracks/thruxton/ThruxtonOfficial.json')
        else:
            print('Not implemented')

        with open(trackFile, 'r') as f:
            self.original_map = json.load(f)
            self.ref_point = self.original_map["ReferencePoint"]

        _out = np.asarray(self.original_map['Outside'])
        _in = np.asarray(self.original_map['Inside'])

        # self.outside_arr = np.asarray(self.original_map['Outside'])[:, :-1]
        # self.inside_arr = np.asarray(self.original_map['Inside'])[:, :-1]
        self.outside_arr = _out if _out.shape[-1] == 2 else _out[:, :-1]
        self.inside_arr = _in if _in.shape[-1] == 2 else _in[:, :-1]
        self.centerline_arr = np.asarray(self.original_map['Centre'])
        
        self._calc_basic_info()

        # Fit a spline to the data - s is the amount of smoothing, tck is the parameters of the             resulting spline
        self.tck_out, uu = si.splprep(self.outside_arr.T, s=0)
        self.tck_in, uu = si.splprep(self.inside_arr.T, s=0)
    
    def _calc_basic_info(self):
        raceline = self.centerline_arr
        self.race_x = raceline[:, 0]
        self.race_y = raceline[:, 1]
        self.raceline_length = self.race_x.shape[0]

        X_diff = np.concatenate([self.race_x[1:] - self.race_x[:-1],
                                 [self.race_x[0] - self.race_x[-1]]])
        Y_diff = np.concatenate([self.race_y[1:] - self.race_y[:-1],
                                 [self.race_y[0] - self.race_y[-1]]])
        self.race_yaw = np.arctan(Y_diff / X_diff)  # (L-1, n)
        
        # arctan returns value in [-pi/2, pi/2)
        self.race_yaw[X_diff < 0] += np.pi
        ## Smooth transition
        self.race_yaw = smooth_yaw(self.race_yaw)
        
        
    def _calc_shortest_distance(self, p, tck, u_init = 0.5):
        # Find the closest point on the spline to our 3d point p
        # We do this by finding a value for the spline parameter u which
        # gives the minimum distance in 3d to p
        # Return distance from 3d point p to a point on the spline at spline parameter u
        def distToP(u):
            s = si.splev(u, tck)
            return ssd.euclidean(p, s)
        output = so.fmin(distToP, u_init, disp = False, full_output = True)
        u_opt = output[0]
        f_opt = output[1]
        p_closest = si.splev(u_opt, tck)
        p_closest = np.array(p_closest).reshape(-1)
        return u_opt.item(), p_closest, f_opt

    def _calc_value(self, p, u_init):
        u_opt, p_out, f_out = self._calc_shortest_distance(p, self.tck_out, u_init = u_init)
        u_opt, p_in, f_in = self._calc_shortest_distance(p, self.tck_in, u_init = u_opt)
        if np.dot(p-p_in, p-p_out) < 0: # inside the track
            v = min(f_out, f_in)
        else: # outside the track
            v = - min(f_out, f_in) 
        return v, u_opt
    
    def get_init_value(self, grid, u_init = 0.5, basis = None): 
        data = np.zeros((grid.pts_each_dim[0], grid.pts_each_dim[1]))
        xx = []
        yy = []
        for idx_x, x in enumerate(grid.vs[0][:, 0, 0, 0]):
            for idx_y, y in enumerate(grid.vs[1][0, :, 0, 0]):
                idx = (idx_x, idx_y)
                if basis is not None:
                    p = basis.dot(np.array([x, y, 1]))
                else:
                    p = np.array([x, y])
                v, _ = self._calc_value(p, u_init)
                data[idx] = v
                xx.append(p[0])
                yy.append(p[1])
        xx = np.array(xx).reshape(data.shape)
        yy = np.array(yy).reshape(data.shape)
        # expand to match the dimension of the grid
        data = np.expand_dims(data, axis = [2, 3])
        data = np.tile(data, (1, 1, grid.pts_each_dim[2], grid.pts_each_dim[3]))
        return data, xx, yy