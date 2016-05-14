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

from medianshape.simplicial import pointgen3d, mesh
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

def surfaces(bbox=[0,0,-10,2*np.pi,2*np.pi,10]):
    '''
    '''
    # Generating point grids for two surfaces
    leftmost= 0
    rightmost  = 2*np.pi
    intersect = np.pi
    step = 0.2
    x = np.arange(leftmost, rightmost, step)
    xb = np.array([leftmost, rightmost])
    x = np.append(x, rightmost)
    y = np.arange(leftmost, rightmost, step)
    yb = np.array([leftmost, rightmost])
    y = np.append(y, rightmost)
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
    s1 = np.concatenate((X.reshape(-1,1), Y.reshape(-1,1), Z1.reshape(-1,1)), axis=1).reshape(-1,3)
    s2 = np.concatenate((X.reshape(-1,1), Y.reshape(-1,1), Z2.reshape(-1,1)), axis=1).reshape(-1,3)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(s1[:,0], s1[:,1], s1[:,2])
    surf = ax.scatter(s2[:,0], s2[:,1], s2[:,2])
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    plt.show()
    leftline1 = np.where(s1[:,0]==leftmost)[0]
    rightline1 = np.where(s1[:, 0] == rightmost)[0]
    backline1 = np.where(s1[:,1]==leftmost)[0]
    frontline1 = np.where(s1[:, 1] == rightmost)[0]
    b1 = np.unique(np.concatenate((leftline1, rightline1, backline1, frontline1), axis=0))
    print b1

    leftline2 = np.where(s2[:,0]==leftmost)[0]
    rightline2 = np.where(s2[:, 0] == rightmost)[0]
    backline2 = np.where(s2[:,1]==leftmost)[0]
    frontline2 = np.where(s2[:, 1] == rightmost)[0]
    b2 = np.unique(np.concatenate((leftline2, rightline2, backline2, frontline2), axis=0))
    intersection = np.where(s1[:,0]== intersect)[0]
    closed_boundary = np.concatenate((leftline1, rightline2), axis=0)
    print b2
    print leftline1
    print rightline1
    print leftline2
    print leftline2
    for p in closed_boundary 
    exit()

    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(backline1[:,0], backline1[:,1], backline1[:,2])
    surf = ax.scatter(frontline1[:,0],frontline1[:,1],frontline1[:,2])
    surf = ax.scatter(leftline1[:,0], leftline1[:,1], leftline1[:,2])
    surf = ax.scatter(rightline1[:,0], rightline1[:,1], rightline1[:,2])
    surf = ax.scatter(backline2[:,0], backline2[:,1], backline2[:,2])
    surf = ax.scatter(frontline2[:,0],frontline2[:,1],frontline2[:,2])
    surf = ax.scatter(leftline2[:,0], leftline2[:,1], leftline2[:,2])
    surf = ax.scatter(rightline2[:,0], rightline2[:,1], rightline2[:,2])
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    plt.show()


    s1_complex = Delaunay(s1[:,:-1])
    s2_complex = Delaunay(s2[:,:-1])
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.plot_trisurf(s1[:,0], s1[:,1], s1[:,2], triangles = s1_complex.simplices, cmap=cm.autumn)
    ax.plot_trisurf(s2[:,0], s2[:,1], s2[:,2], triangles = s2_complex.simplices, cmap=cm.winter)
    plt.show()
    exit()

    s_points = np.vstack((s1, s2))
    s1_triangles = s1_complex.simplices
    s2_triangles = s2_complex.simplices + len(s1_complex.points)
    surfaces = np.vstack((s1_triangles, s2_triangles))

    # Plotting the surfaces
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.plot_trisurf(s_points[:,0], s_points[:,1], s_points[:,2], triangles = s1_triangles, cmap=cm.autumn)
    ax.plot_trisurf(s_points[:,0], s_points[:,1], s_points[:,2], triangles = s2_triangles, cmap=cm.winter)
    plt.show()

    '''
    # Plotting the surfaces
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.set_zlim(-10, 10)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles = surfaces, cmap=cm.autumn)
    plt.show()
    '''
    # Grid Points sampled from boundary box
    corners = boundary_points(bbox) 
    box = np.array(bbox).reshape(2, -1)
    points = np.mgrid[tuple(slice(min, max+1, 1) for min, max in box.T)]
    points = points.reshape(3, -1).T

    # Points of PLC.
    # Uncomment the line below to include sample poinds in bbox along with the surface points  
    #points  = np.concatenate((s_points, points), axis=0).reshape(-1,3)
    points  = np.concatenate((s_points, corners), axis=0).reshape(-1,3)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(points[:,0], points[:,1], points[:,2])
    ax.set_zlim(-12, 12)
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    plt.show()

    # Saving the current PLC to a file in .smesh format
    with open('surf.smesh', 'w') as f:
        f.write('%d %d %d %d\n'%(len(points), 3, 0,0))
        for i, p in enumerate(points):
            f.write("%d %f %f %f\n"%(i, p[0], p[1], p[2]))
        f.write('%d %d\n'%(len(surfaces), 0))
        for t in surfaces:
            f.write("%d %d %d %d\n"%(3, t[0], t[1], t[2]))
        f.write("%d\n"%0)
        f.write("%d\n"%0)
    np.savetxt("surface.face", surfaces)
    np.savetxt("points.node", points)

    # Build the mesh using Tetgen
    mesh_info = MeshInfo()
    mesh_info.set_points(points.tolist())
    mesh_info.set_facets((surfaces.reshape(-1, 3)).tolist())
    print len(mesh_info.facets)
    opts = Options("YVfeq", verbose=True, nobisect=True, facesout=True, edgesout=True, docheck=True) # Overriding 'pq' with no options or flags
    #opts = Options("", verbose=True, nobisect=True, facesout=True, edgesout=True, docheck=True, insertaddpoints=True) # Overriding 'pq' with no options or flags
    mesh = build(mesh_info, options=opts, volume_constraints=True, max_volume=1)

    # Plot the mesh
    ax = plt.gca(projection='3d')
    fig = ax.figure
    m1 = np.amin(bbox[0:3])
    m2 = np.amax(bbox[3:])
    ax.set_xlim([m1, m2])
    ax.set_ylim([m1, m2])
    ax.set_zlim([m1, m2])
    ax.set_aspect('equal')
    axes_simpplot3d(ax, np.array(list(mesh.points)), np.array(list(mesh.elements)))
    plt.show()
    # Write it as a file in vtk format, so you can use Paraview to see it.
    mesh.write_vtk("test.vtk")
    input_surfaces = np.zeros((2,len(mesh.faces)))
    inputs = list (s1_triangles, s2_triangles)
    for i, s in enumerate(inputs):
        for j, t in enumerate(np.array(mesh.faces)):
            if np.all(s==t):
                input_surfaces[i,j] = 1

    lambdas = [0.001]
    mus = [0.00001]
    mesh1 = Mesh3D()
    mesh1.simplices = np.array(mesh.elements)
    mesh1.triangles = np.array(mesh.faces)
    mesh1.edges = np.array(mesh.edges)
    mesh1.points = np.array(mesh.points)
    return mesh1, mesh1.simplices, mesh1.triangles, input_surfaces, lambdas, mus

    
if __name__ == "__main__":
    surfaces()
