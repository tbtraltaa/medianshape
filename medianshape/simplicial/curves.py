# encoding: utf-8
'''
Utils
=====
'''

from __future__ import division

import sys
import importlib

import numpy as np
import math
import matplotlib.pyplot as plt

def curve1(x):
    '''
    :math:`f(x) = 10e^{\\frac{x}{180}}\sin(\\frac{1}{100}\pi x) +25`
    '''
    return 10*math.sin(1/100* math.pi * x) * np.exp(1/180 * x)  + 25

def curve2(x):
    '''
    :math:`f(x) = -10e^{\\frac{x}{200}}\sin(\\frac{1}{50}\pi x) +25`
    '''
    return -10*math.sin(1/50* math.pi * x) * np.exp(1/200 * x) + 25

def curve3(x):
    '''
    :math:`f(x) = 30e^{\\frac{-x}{150}}\sin(\\frac{1}{50}\pi x) +20`
    '''
    return 30*math.sin(1/50* math.pi * x) * np.exp(-1/150 * x) + 25
    
if __name__ == "__main__":
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
        #Ys4.append(curve4(x))
        #Ys5.append(curve5(x))
        #Ys6.append(curve6(x))
    #plt.ylim(0, 60)
    #plt.xlim(0, 210)
    plt.figure()
    #plt.plot(Xs, Ys, '-', label="1")
    #plt.plot(Xs, Ys2, '-', label="2")
    #plt.plot(Xs, Ys3, '-', label="3")
    plt.plot(Xs, Ys, '-', label="1", color='r')
    plt.plot(Xs, Ys2, '-', color='b', label="2")
    plt.plot(Xs, Ys3, '-', color='y', label="3")
    plt.show()
