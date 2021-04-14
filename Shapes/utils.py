import math
import numpy as np

def toLocal(XY, origin, yaw):
    # (x, y) -> (x', y')
    R = np.array([[np.cos(yaw), np.sin(yaw)], 
                  [-np.sin(yaw), np.cos(yaw)]])
    return (XY-origin).dot(R.T)

def toGlobal(XY, origin, yaw):
    # (x', y') -> (x, y)
    R_inv = np.array([[np.cos(yaw), -np.sin(yaw)], 
                  [np.sin(yaw), np.cos(yaw)]])
    return XY.dot(R_inv.T) + origin

def fitRectangle(pts):
    x_min = np.infty
    x_max = -np.infty
    y_min = np.infty
    y_max = -np.infty

    for x, y in pts:
        x_min = min(x_min, x)
        x_max = max(x_max, x)
        y_min = min(y_min, y)
        y_max = max(y_max, y)

    return x_min, x_max, y_min, y_max

def smooth_yaw(yaw):
    for i in range(len(yaw) - 1):
        dyaw = yaw[i + 1] - yaw[i]

        while dyaw >= math.pi / 2.0:
            yaw[i + 1] -= math.pi * 2.0
            dyaw = yaw[i + 1] - yaw[i]

        while dyaw <= -math.pi / 2.0:
            yaw[i + 1] += math.pi * 2.0
            dyaw = yaw[i + 1] - yaw[i]
    return yaw