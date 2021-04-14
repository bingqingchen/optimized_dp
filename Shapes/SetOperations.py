import numpy as np

def setUnion(data1, data2):
    assert data1.shape == data2.shape
    return np.minimum(data1, data2)

def setIntersection(data1, data2):
    assert data1.shape == data2.shape
    return np.maximum(data1, data2)

def setComplement(data):
    return -data

def setDifference(data1, data2):
    # S1\S2
    assert data1.shape == data2.shape
    return np.maximum(data1, -data2)


