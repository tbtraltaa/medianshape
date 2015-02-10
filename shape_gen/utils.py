# encoding: utf-8

from __future__ import division

import sys
import importlib
import math

import numpy as np
import matplotlib.pyplot as plt

def get_combination(n, vals=[1,-1]):
    if n > 1:
        comb = np.array([])
        successors = get_combination(n-1, vals)
        for val in vals:
            for s in successors:
                row = np.append(val, s)
                if comb.size == 0:
                    comb = row
                else:
                    comb = np.vstack((comb,row))
        return comb
    elif n==1:
        return vals

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
    return 5*np.abs(math.sin(math.pi*x)) + 10
def half_sin1pi(x):
    return 2*np.abs(math.sin(math.pi*x))/2 + 5
def small_sin1pi(x):
    return np.abs(math.sin(math.pi*x))/4
def func2(x):
    return 1/2*(1+math.sin(2*math.pi*x))

#def curve1(x, boundary_box):
#    no_of_intervals = 8
#    heights = [5, 10, 5, 5, 10, 5, 10, 5]
#    offsets = [10, 20, 10, 10, 20, 10, 20, 10]
#    phases = [1, 2, 1, 1, 2, 1, 2, 1]
#    numbers = np.linspace(boundary_box[0], boundary_box[2], no_of_intervals) 
#    intervals = [(a, numbers[i+1]) for i, a in enumerate(numbers) if i < len(numbers) - 1]
#    interval_number = [i for i, interval in enumerate(intervals) if x >= interval[0] and x < interval[1]]
#    j = interval_number[0]
#    return heights[j]*math.sin(phases[j]*math.pi*x) + offsets[j]
# green
def curve1(x):
    return 10*math.sin(1/100* math.pi * x) * np.exp(1/180 * x)  + 25
# red
def curve2(x):
    return -10*math.sin(1/50* math.pi * x) * np.exp(1/200 * x) + 25
# violate
def curve3(x):
    return 30*math.sin(1/50* math.pi * x) * np.exp(-1/150 * x) + 20

def curve4(x):
    return -40*math.sin(1/200* math.pi * x) + 45

def curve5(x):
    return 10*math.sin(1/30* math.pi * x) * np.exp(-1/200 * x) + 10
def curve6(x):
    return 10*math.sin(1/30* math.pi * x) * np.exp(-1/200 * x) + 40

    
if __name__ == "__main__":
    print  get_combination(3)
    Xs = np.linspace(0, 200, 20)
    Ys = list()
    Ys2 = list()
    Ys3 = list()
    Ys4 = list()
    Ys5 = list()
    Ys6 = list()
    for x in Xs:
        Ys.append(curve1(x))
        Ys2.append(curve2(x))
        Ys3.append(curve3(x))
        Ys4.append(curve4(x))
        Ys5.append(curve5(x))
        Ys6.append(curve6(x))
    #plt.ylim(0, 60)
    #plt.xlim(0, 210)
    plt.figure()
    plt.plot(Xs, Ys, '-', label="1")
    plt.plot(Xs, Ys2, '-', label="2")
    plt.plot(Xs, Ys3, '-', label="3")
    plt.plot(Xs, Ys4, '-', label="4")
    plt.plot(Xs, Ys5, '-', label="5")
    plt.plot(Xs, Ys6, '-', label="6")
    #plt.show()