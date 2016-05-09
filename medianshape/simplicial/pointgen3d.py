# encoding: utf-8
'''
Point generation 3D
===================
'''

from __future__ import division

import numpy as np

def curve1(bbox):
    '''
    Hi
    '''
    x = np.linspace(bbox[0], bbox[3], 5).reshape(-1, 1)
    points = np.tile(x,(1,3))
    return points

def curve2(bbox):
    '''
    Hi
    '''
    x = np.linspace(bbox[0], bbox[3], 5).reshape(-1, 1)
    z = x**10
    points = np.tile(x,(1,2))
    points = np.hstack((points,z))
    return points

def sphere_arc(bbox, theta, n):
    '''
    Hi
    '''
    r = np.abs(bbox[3] - bbox[0])*1.0/2  
    center = [(bbox[0]+bbox[3])*1.0/2, (bbox[1] + bbox[4])*1.0/2, (bbox[2]+bbox[5])*1.0/2]
    verticle_angle = np.linspace(0, np.pi, n)
    x = r*np.cos(theta)*np.sin(verticle_angle) + center[0]
    y = r*np.sin(theta)*np.sin(verticle_angle) + center[1]
    z = r*np.cos(verticle_angle) + center[2]
    return np.hstack((x.reshape(-1,1), y.reshape(-1,1), z.reshape(-1,1)))

def horizontal_circle(r, n=10, center=[0,0,0], theta=0):
    x = center[0] + r*np.cos(np.linspace(0, 2*np.pi, n)).reshape(-1,1)
    y = center[1] + r*np.sin(np.linspace(0, 2*np.pi, n)).reshape(-1,1)
    z = center[2] + np.zeros(n).reshape(-1,1)
    return np.concatenate((x,y,z), axis=1)

def vertical_circle_xz(r, n=5, center=[0,0,0], theta=0):
    x = center[0] + r*np.cos(np.linspace(0, 2*np.pi, n)).reshape(-1,1)
    y = center[1] + np.zeros(n).reshape(-1,1)
    z = center[2] + r*np.sin(np.linspace(0, 2*np.pi, n)).reshape(-1,1)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    if theta != 0:
        x1 = x*cos_theta - y*sin_theta
        y1 = y*cos_theta + x*sin_theta
        x = x1
        y = y1
    return np.concatenate((x,y,z), axis=1)

if __name__ == "__main__":
    curve1([0,0,0, 10,10,10])
