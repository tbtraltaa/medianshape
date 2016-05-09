# encoding: utf-8
'''
2D Median surface embedded in 3D
--------------------------------
'''
from __future__ import absolute_import
import importlib
import os

import numpy as np
from scipy.spatial import Delaunay

from medianshape.simplicial import pointgen3d
from medianshape.simplicial.meshgen import meshgen3d, get_mesh_surface
import medianshape.experiment.inout as inout

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import cm

import medianshape.viz.plot3d as plot3d
from distmesh.plotting import axes_simpplot3d
from meshpy.tet import MeshInfo, Options, build
from medianshape.simplicial.utils import boundary_points

def surfaces(bbox=[0,0,-10,10,10,10]):
    '''
    '''
    # Generate points in a cube given as a bbox
    '''
    bbox = np.array(bbox).reshape(2, -1)
    dim = bbox.shape[1]
    points = np.mgrid[tuple(slice(min, max+l, l) for min, max in bbox.T)]
    points = points.reshape(dim, -1).T
    '''
    # Generating point grids for two surfaces
    x = np.arange(0, 3*np.pi, 0.3)
    x = np.append(x, 3*np.pi)
    y = np.arange(0, 3*np.pi, 0.3)
    y = np.append(y, 3*np.pi)
    X, Y = np.meshgrid(x, y, sparse=False)
    # Z coordinate of surface1
    Z1 = 7*np.sin(X)
    # Z coordinate of surface2
    Z2 = -7*np.sin(X)
    # Illustrating the surfaces
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, cmap=cm.autumn,
                           linewidth=0, antialiased=False)
    surf = ax.plot_surface(X, Y, Z2, rstride=1, cstride=1, cmap=cm.winter,
                           linewidth=0, antialiased=False)
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    plt.show()

    # Triangulating the surfaces
    s1 = np.concatenate((X.reshape(-1,1), Y.reshape(-1,1), Z1.reshape(-1,1)), axis=1)
    s2 = np.concatenate((X.reshape(-1,1), Y.reshape(-1,1), Z2.reshape(-1,1)), axis=1)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(s1[:,0], s1[:,1], s1[:,2])
    surf = ax.scatter(s2[:,0], s2[:,1], s2[:,2])
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    plt.show()

    s1_complex = Delaunay(s1)
    s2_complex = Delaunay(s2)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.plot_trisurf(s1[:,0], s1[:,1], s1[:,2], triangles = s1_complex.simplices, cmap=cm.autumn)
    ax.plot_trisurf(s2[:,0], s2[:,1], s2[:,2], triangles = s2_complex.simplices, cmap=cm.winter)
    plt.show()

    corners = boundary_points(bbox) 
    box = np.array(bbox).reshape(2, -1)
    points = np.mgrid[tuple(slice(min, max+1, 1) for min, max in box.T)]
    points = points.reshape(3, -1).T
    points  = np.concatenate((corners, points, s1,s2), axis=0).reshape(-1,3)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(points[:,0], points[:,1], points[:,2])
    surf = ax.scatter(points[:,0], points[:,1], points[:,2])
    ax.set_zlim(-12, 12)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    plt.show()
    surface = np.hstack((s1_complex.simplices, s2_complex.simplices))
    print "hi"
    #print len(mesh_info.points)
    mesh_info = MeshInfo()
    mesh_info.set_points(points)
    mesh_info.set_facets(surface)
    #opts = Options(switches="pY", p=points, Y=surface) # Overriding 'pq' with no options or flags
    opts = Options(switches="") # Overriding 'pq' with no options or flags
    mesh = build(mesh_info, options=opts, volume_constraints=True, max_volume=3)
    print "finished"
    ax = plt.gca(projection='3d')
    fig = ax.figure
    m1 = np.amin(bbox[0:3])
    m2 = np.amax(bbox[3:])
    ax.set_xlim([m1, m2])
    ax.set_ylim([m1, m2])
    ax.set_zlim([m1, m2])
    ax.set_aspect('equal')
    #axes_simpplot3d(ax, mesh.points, mesh.simplices, mesh.points[:,1] > 0)
    print len(mesh.elements)
    for e in mesh.elements:
        print e
        
    axes_simpplot3d(ax, np.array(list(mesh.points)), np.array(list(mesh.elements)))
    plt.show()
    for f in surface:
        if f 
    
if __name__ == "__main__":
    surfaces()
