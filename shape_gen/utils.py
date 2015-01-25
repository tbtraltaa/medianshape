# encoding: utf-8

from __future__ import division

import sys
import importlib
import math

import numpy as np

def vectorize(func_str, X):
    func_points = []
    if func_str.find(".") != -1:
        mod_name, func_name = func_str.rsplit('.', 1)
        mod = importlib.import_module(mod_name)
        func = getattr(mod, func_name, "")
        if func:
            vec_func = np.vectorize(func)    
    else:
        this_module = sys.modules[__name__]
        func = getattr(this_module, func_str, "")
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
