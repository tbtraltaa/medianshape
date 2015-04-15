# encoding: utf-8

from __future__ import division, absolute_import

import importlib
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, cdist


def func_sampling(func_str=None, boundary_box, sample_size, mesh=None):
        sample_X = []
        if mesh:
            ordered_X = np.sort(np.unique(mesh.points[:,0]))
            sample_size = int(math.sqrt(len(mesh.points)))
            sample_step = ordered_X.size//sample_size
            if sample_step == 0:
                sample_step = 1
            for i in range(0, ordered_X.size, self.sample_step):
                sample_X.append(ordered_X[i])
        else:
            sample_X = np.logspace(boundary_box[0], boundary_box[2], sample_size)
        func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), \
        vectorize(func_str, sample_X).reshape(len(sample_X),1)))
        return func_points

def vectorize(func_str, X):
    func_points = []
    if func_str.find(".") != -1:
        mod_name, func_name = func_str.rsplit('.', 1)
        mod = importlib.import_module(mod_name)
        func = getattr(mod, func_name)
        vec_func = np.vectorize(func)    
    else:
        func = getattr(FunctionApprox2d, func_str)
        vec_func = np.vectorize(func)    
    func_values = vec_func(X)
    return func_values

def x(x):
    return x
def x2(x):
    return x**2
def x5(x):
    return x**5
def func1(x):
    return 2/math.pi*math.acos(x)
def sin2pi(x):
    return np.abs(math.sin(2*math.pi*x))
def sin1pi(x):
    return np.abs(math.sin(math.pi*x))
def half_sin1pi(x):
    return np.abs(math.sin(math.pi*x))/2
def small_sin1pi(x):
    return np.abs(math.sin(math.pi*x))/4
def func2(x):
    return 1/2*(1+math.sin(2*math.pi*x))
