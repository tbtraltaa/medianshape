# encoding: utf-8

from __future__ import division

import numpy as np

def curve1(bbox):
    x = np.linspace(bbox[0], bbox[3], 5).reshape(-1, 1)
    points = np.tile(x,(1,3))
    return points

def curve2(bbox):
    x = np.linspace(bbox[0], bbox[3], 5).reshape(-1, 1)
    z = x**10
    points = np.tile(x,(1,2))
    points = np.hstack((points,z))
    return points

def sphere_arc(bbox, theta, n):
    r = np.abs(bbox[3] - bbox[0])*1.0/2  
    center = [(bbox[0]+bbox[3])*1.0/2, (bbox[1] + bbox[4])*1.0/2, (bbox[2]+bbox[5])*1.0/2]
    verticle_angle = np.linspace(0, np.pi, n)
    x = r*np.cos(theta)*np.sin(verticle_angle) + center[0]
    y = r*np.sin(theta)*np.sin(verticle_angle) + center[1]
    z = r*np.cos(verticle_angle) + center[2]
    return np.hstack((x.reshape(-1,1), y.reshape(-1,1), z.reshape(-1,1)))

def sphere_equator(bbox, n):
    r = np.abs(bbox[3] - bbox[0])*1.0/2  
    center = [(bbox[0]+bbox[3])*1.0/2, (bbox[1] + bbox[4])*1.0/2, (bbox[2]+bbox[5])*1.0/2]
    theta = np.linspace(0, np.pi, n)
    x = r*np.cos(theta) + center[0]
    y = r*np.sin(theta) + center[1]
    z = np.tile(center[2],(n,1))
    return np.hstack((x.reshape(-1,1), y.reshape(-1,1), z.reshape(-1,1)))



if __name__ == "__main__":
    curve1([0,0,0, 10,10,10])
