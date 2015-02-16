# encoding: utf-8

from __future__ import division

import sys
import importlib
import math

import numpy as np

import utils

# Samples function with the boundary box
def sample_function(func_str, boundary_box, sample_size):
        sample_X = np.linspace(boundary_box[0], boundary_box[2], sample_size)
        func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), \
        utils.vectorize(func_str, sample_X).reshape(len(sample_X),1)))
        return func_points

# Samples function based on the mesh given
def sample_function_mesh(func_str, mesh, sample_size=None):
        sample_X = []
        ordered_X = np.sort(np.unique(mesh.points[:,0]))
        if not sample_size:
            sample_size = int(math.sqrt(len(ordered_X)))
        sample_step = ordered_X.size//sample_size
        if sample_step == 0:
            sample_step = 1
        for i in range(0, ordered_X.size, sample_step):
            sample_X.append(ordered_X[i])
        func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), \
        utils.vectorize(func_str, sample_X).reshape(len(sample_X),1)))
        return func_points

def sample_ellipse(a, b, sample_size=None):
    return np.array([(a * math.cos(theta) +0.5 , b * math.sin(theta)+0.5)\
            for theta in (math.pi*2 * i/sample_size for i in range(sample_size))])

if __name__ == "__main__":
    print sample_function_on_boundary("x5", (0,0,1,1), 20);
