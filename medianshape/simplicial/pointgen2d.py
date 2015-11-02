# encoding: utf-8
'''
Point generation 2D
===================

'''

from __future__ import division

import numpy as np
import math

from medianshape.simplicial.utils import vectorize

# Samples function with the boundary box
def sample_function(func_str, boundary_box, sample_size):
    '''
    Hi
    '''
    sample_X = np.linspace(boundary_box[0], boundary_box[2], sample_size)
    func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), \
    vectorize(func_str, sample_X).reshape(len(sample_X),1)))
    return func_points

# Samples function based on the mesh given
def sample_function_mesh(mesh,  func_str, sample_size=None):
    '''
    Hi
    '''
    sample_X = []
    ordered_X = np.sort(np.unique(mesh.points[:,0]))
    if not sample_size:
        sample_size = int(math.sqrt(len(mesh.points)))
    sample_step = ordered_X.size//sample_size
    if sample_step == 0:
        sample_step = 1
    for i in range(0, ordered_X.size, sample_step):
        sample_X.append(ordered_X[i])
    func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), \
    vectorize(func_str, sample_X).reshape(len(sample_X),1)))
    return func_points

def sample_ellipse(a, b, sample_size=None):
    '''
    Hi
    '''
    return np.array([(a * math.cos(theta) +0.5 , b * math.sin(theta)+0.5)\
            for theta in (math.pi*2 * i/sample_size for i in range(sample_size))])

def twisted_curve1():
    '''
    Hi
    '''
    x = np.linspace(0,30, 5)
    x = np.append(x, np.linspace(30, 0, 5)[1:])
    x = np.append(x, np.linspace(0, 35, 5)[1:])
    x = np.append(x, np.linspace(35, 30, 5)[1:])
    y = np.linspace(0,15, 5)
    y = np.append(y, np.linspace(15, 20, 5)[1:])
    y = np.append(y, np.linspace(20, 25, 5)[1:])
    y = np.append(y, np.linspace(25, 0, 5)[1:])
    x = x.reshape((x.shape[0], 1))
    y = y.reshape((y.shape[0], 1))
    points = np.hstack((x,y))
    return points

def twisted_curve2():
    '''
    Hi
    '''
    x = np.linspace(0,15, 5)
    x = np.append(x, np.linspace(15,20, 5)[1:])
    x = np.append(x, np.linspace(20, 30, 5)[1:])
    y = np.linspace(0,30, 5)
    y = np.append(y, np.linspace(30, 5, 5)[1:])
    y = np.append(y, np.linspace(5, 0, 5)[1:])
    x = x.reshape((x.shape[0], 1))
    y = y.reshape((y.shape[0], 1))
    points = np.hstack((x,y))
    return points

if __name__ == "__main__":
    #print sample_function_on_boundary("x5", (0,0,1,1), 20);
    twisted_curve1()
    twisted_curve2()
