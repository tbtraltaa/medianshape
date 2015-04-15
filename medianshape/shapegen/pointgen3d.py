# encoding: utf-8

from __future__ import division

import sys
import importlib
import math

import numpy as np
import matplotlib.pyplot as plt

import utils

def curve1(bbox):
    x = np.linspace(bbox[0], bbox[3], 5).reshape(-1, 1)
    points = np.tile(x,(1,3))
    return points

def curve2(bbox):
    x = np.linspace(bbox[0], bbox[3], 5).reshape(-1, 1)
    z = x**5
    points = np.tile(x,(1,2))
    points = np.hstack((points,z))
    return points


if __name__ == "__main__":
    curve1([0,0,0, 10,10,10])
