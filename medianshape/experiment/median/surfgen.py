# encoding: utf-8
'''
2D Median surface embedded in 3D
--------------------------------
'''
from __future__ import absolute_import
import importlib
import os

import numpy as np

from medianshape.simplicial import pointgen3d, mesh, utils
from medianshape.simplicial.meshgen import meshgen2d

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from medianshape.viz import plot2d, plot3d
from distmesh.plotting import axes_simpplot3d
from medianshape.simplicial.utils import boundary_points

def func(x, y, sign=1): 
    return np.sin(np.pi*x)*np.cos(np.pi*y)

def sample_surf(scale, step=0.2):
    '''
        Return a tuple X, Y, Z with a test surface.
    '''
    x = y = np.arange(-4.0, 4.0, step)
    X, Y = np.meshgrid(x, y)
    from matplotlib.mlab import  bivariate_normal
    '''
    Z1 = bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
    Z2 = bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
    #Z3 = bivariate_normal(X, Y, 1, 1, -2, -2)
    Z = Z2 - Z1
    '''
    # Ups
    ZU1 = bivariate_normal(X,Y, 1.5, 1, 0,-2)
    ZU2 = bivariate_normal(X, Y, 1.5, 1.5, 4, 1)
    ZU3 = bivariate_normal(X, Y, 1, 1, -4, 1)
    #ZU4 = bivariate_normal(X, Y, 1.5, 1.5, -4, -4)
    #ZU5 = bivariate_normal(X, Y, 1, 1, 4, -4)
    ZU4 = bivariate_normal(X, Y, 4, 0.5, 0, -4)
    # Downs
    ZD1 = bivariate_normal(X, Y, 1.5, 1, 0, 1)
    ZD2 = bivariate_normal(X, Y, 1.5, 1.5, -4, -2)
    ZD3 = bivariate_normal(X, Y, 1, 1, 4, -2)
    ZD4 = bivariate_normal(X, Y, 4, 1, 0, 4)
    Z1 = ZU1 + ZU2 + ZU3 - ZD1 - ZD2 - ZD3 - ZD4
    Zmax1 = np.abs(np.amax(Z1))
    Z1 = Z1/Zmax1 * scale[2]
    
    # Visualization
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, cmap=cm.winter,
                           linewidth=0, antialiased=False)
    plt.show()

    # Ups
    ZU1 = bivariate_normal(X,Y, 2, 1, 1,1)
    ZU2 = bivariate_normal(X, Y, 3, 1, -2, 4)
    ZU3 = bivariate_normal(X, Y, 1.5, 1.5, -2, -2)
    #ZU4 = bivariate_normal(X, Y, 1.5, 1.5, -4, -4)
    #ZU5 = bivariate_normal(X, Y, 1, 1, 4, -4)
    ZU4 = bivariate_normal(X, Y, 2, 2, 3, -4)
    # Downs
    ZD1 = bivariate_normal(X, Y, 1, 2, 4, 2)
    ZD2 = bivariate_normal(X, Y, 1.5, 1.5, -2, 2)
    ZD3 = bivariate_normal(X, Y, 1.5, 1.5, 1, -2)
    ZD4 = bivariate_normal(X, Y, 4, 1, 0, -4)
    Z2 = ZU1 + ZU2 + ZU3 - ZD1 - ZD2 - ZD3 - ZD4
    Zmax2 = np.abs(np.amax(Z2))
    Z2 = Z2/Zmax2 * scale[2]

    X = X * scale[0]/4.0
    Y = Y * scale[1]/4.0

    # Visualization
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.plot_surface(X, Y, Z2, rstride=1, cstride=1, cmap=cm.winter,
                           linewidth=0, antialiased=False)
    plt.show()
    return X, Y, Z1, Z2


def interpolate_surf(points, values, ipoints, method = "nearest"):
    from scipy.interpolate import griddata
    return griddata(points, values, ipoints, method= method)

def surfgen_shared_boundary(bbox=[-10,-10,-10, 10,10,10], l=3):
    '''
    '''
    # Generating point grids for two surfaces
    xmin = bbox[0]
    xmax  = bbox[3]
    ymin = bbox[1]
    ymax = bbox[4]
    zmin = bbox[2]
    zmax = bbox[5]
    Xmin, Ymin, Zmin, Xmax, Ymax, Zmax = np.array(bbox)*0.8
    X, Y, Z1, Z2 = sample_surf([Xmax, Ymax, zmax*0.3], step=0.8)
    Z1 = Z1 + zmax*0.4
    Z2 = Z2 - zmax*0.4
    #Symmertic surfs
    #Z2 = -Z1 - zmax*0.4
    '''
    # Plotting the two surfaces
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(X, Y, Z1.reshape(-1,1), color='b')
    surf = ax.scatter(X, Y, Z2.reshape(-1,1), color='r')
    plt.show()
    '''
    mesh = meshgen2d([Xmin, Ymin, Xmax, Ymax], l, include_corners=True)
    sample_points = np.hstack((X.reshape(-1,1), Y.reshape(-1,1)))

    # Interpolating the surface mesh into two different surfaces
    # similar to the the sample surfaces generated before
    Z1 = interpolate_surf(sample_points, Z1.reshape(-1,1), mesh.points)
    Z2 = interpolate_surf(sample_points, Z2.reshape(-1,1), mesh.points)
    
    # Integrating two surfaces
    points1 = np.hstack((mesh.points, Z1))
    print points1.shape
    points2 = np.hstack((mesh.points, Z2))
    print points2.shape
    corners = utils.boundary_points(bbox)
    midcorners = utils.mid_corners(bbox)
    offset1 = len(corners) +len(midcorners) + 1
    offset2 = len(corners) + len(midcorners) + len(points1) + 1
    points = np.concatenate((corners, midcorners, points1, points2), axis=0)
    print points.shape

    triangles1 = mesh.simplices + offset1
    triangles2 = mesh.simplices + offset2
   
    # Adding the indices of the points as the last column of the coordainate list
    Xmin_s1 = np.argwhere(points1[:,0]==Xmin)
    Xmin_s1_points = np.hstack((points1[Xmin_s1.reshape(-1,)], Xmin_s1))
    # Sorting the indices such that the points are in increasing order of its y-component
    Xmin_s1 = (Xmin_s1_points[:,3][np.argsort(Xmin_s1_points[:,1])] + offset1).astype(int)

    Xmin_s2 = np.argwhere(points2[:,0]==Xmin)
    Xmin_s2_points = np.hstack((points2[Xmin_s2.reshape(-1,)], Xmin_s2))
    Xmin_s2 = (Xmin_s2_points[:,3][np.argsort(Xmin_s2_points[:,1])] + offset2).astype(int)

    Xmax_s1 = np.argwhere(points1[:,0]==Xmax)
    Xmax_s1_points = np.hstack((points1[Xmax_s1.reshape(-1,)], Xmax_s1))
    Xmax_s1 = (Xmax_s1_points[:,3][np.argsort(Xmax_s1_points[:,1])] + offset1).astype(int)

    Xmax_s2 = np.argwhere(points2[:,0]==Xmax)
    Xmax_s2_points = np.hstack((points2[Xmax_s2.reshape(-1,)], Xmax_s2))
    Xmax_s2 = (Xmax_s2_points[:,3][np.argsort(Xmax_s2_points[:,1])] + offset2).astype(int)

    Ymin_s1 = np.argwhere(points1[:,1]==Ymin)
    Ymin_s1_points = np.hstack((points1[Ymin_s1.reshape(-1,)], Ymin_s1))
    Ymin_s1 = (Ymin_s1_points[:,3][np.argsort(Ymin_s1_points[:,0])] + offset1).astype(int)

    Ymin_s2 = np.argwhere(points2[:,1]==Ymin)
    Ymin_s2_points = np.hstack((points2[Ymin_s2.reshape(-1,)], Ymin_s2))
    Ymin_s2 = (Ymin_s2_points[:,3][np.argsort(Ymin_s2_points[:,0])] + offset2).astype(int)

    Ymax_s1 = np.argwhere(points1[:,1]==Ymax)
    Ymax_s1_points = np.hstack((points1[Ymax_s1.reshape(-1,)], Ymax_s1))
    Ymax_s1 = (Ymax_s1_points[:,3][np.argsort(Ymax_s1_points[:,0])] + offset1).astype(int)

    Ymax_s2 = np.argwhere(points2[:,1]==Ymax)
    Ymax_s2_points = np.hstack((points2[Ymax_s2.reshape(-1,)], Ymax_s2))
    Ymax_s2 = (Ymax_s2_points[:,3][np.argsort(Ymax_s2_points[:,0])] + offset2).astype(int)

    for i in range(len(Xmin_s1)-1):
        triangles1 = np.vstack((triangles1, [9, Xmin_s1[i], Xmin_s1[i+1]]))
    triangles1 = np.vstack((triangles1, [9, Xmin_s1[-1], 12]))
    for i in range(len(Xmin_s2)-1):
        triangles2 = np.vstack((triangles2, [9, Xmin_s2[i], Xmin_s2[i+1]]))
    triangles2 = np.vstack((triangles2, [9, Xmin_s2[-1], 12]))
    for i in range(len(Xmax_s1)-1):
        triangles1 = np.vstack((triangles1, [10, Xmax_s1[i], Xmax_s1[i+1]]))
    triangles1 = np.vstack((triangles1, [10, Xmax_s1[-1], 11]))
    for i in range(len(Xmax_s2)-1):
        triangles2 = np.vstack((triangles2, [10, Xmax_s2[i], Xmax_s2[i+1]]))
    triangles2 = np.vstack((triangles2, [10, Xmax_s2[-1], 11]))

    for i in range(len(Ymin_s1)-1):
        triangles1 = np.vstack((triangles1, [9, Ymin_s1[i], Ymin_s1[i+1]]))
    triangles1 = np.vstack((triangles1, [9, Ymin_s1[-1], 10]))
    for i in range(len(Ymin_s2)-1):
        triangles2 = np.vstack((triangles2, [9, Ymin_s2[i], Ymin_s2[i+1]]))
    triangles2 = np.vstack((triangles2, [9, Ymin_s2[-1], 10]))
    
    for i in range(len(Ymax_s1)-1):
        triangles1 = np.vstack((triangles1, [12, Ymax_s1[i], Ymax_s1[i+1]]))
    triangles1 = np.vstack((triangles1, [12, Ymax_s1[-1], 11]))
    for i in range(len(Ymax_s2)-1):
        triangles2 = np.vstack((triangles2, [12, Ymax_s2[i], Ymax_s2[i+1]]))
    triangles2 = np.vstack((triangles2, [12, Ymax_s2[-1], 11]))

    triangles = np.vstack((triangles1, triangles2)) 

    print triangles1.shape
    print triangles2.shape
    print triangles.shape
    
    # Preparing PLC and save it to .poly file for tetgen
    with open( os.environ['HOME'] +'/mediansurf.poly', 'w') as f:
        f.write("#Part 1 - the node list\n")
        f.write("#%d nodes in 3d, no attributes, no boundary marker\n"%points.shape[0])
        f.write('%d %d %d %d\n'%(points.shape[0], 3, 0,0))
        for i, p in enumerate(points):
            f.write("%d %f %f %f\n"%(i+1, p[0], p[1], p[2]))
        # Each 4 sides has 3 polygons
        # Top and bottom
        # Each triangle of the two surfaces are facets
        fn = 6 + len(triangles)
        f.write("#Part 2 - the facet list.\n")
        f.write("#%d facets with boundary markers\n"%fn)
        f.write('%d %d\n'%(fn, 1))
        f.write("#Boundary facet list.\n")
        f.write("%d %d %d\n"%(1, 0, 1))
        f.write("4 1 2 3 4\n")
        f.write("%d %d %d\n"%(1, 0, 1))
        f.write("4 5 6 7 8\n")
        #xmin side
        f.write("2 0 1\n")
        f.write("4 1 4 8 5\n")
        f.write("2 9 12\n")
        #ymin side
        f.write("2 0 1\n")
        f.write("4 1 2 6 5\n")
        f.write("2 9 10\n")
        #xmax side
        f.write("2 0 1\n")
        f.write("4 2 3 7 6\n")
        f.write("2 10 11\n")
        #ymax side
        f.write("2 0 1\n")
        f.write("4 3 4 8 7\n")
        f.write("2 11 12\n")

        f.write("#Facet list of surface1.\n")
        for t in triangles1:
            f.write("%d %d %d\n"%(1, 0, -1))
            f.write("%d %d %d %d\n"%(3, t[0], t[1], t[2]))
        f.write("#Facet list of surface2.\n")
        for t in triangles2:
            f.write("%d %d %d\n"%(1, 0, -2))
            f.write("%d %d %d %d\n"%(3, t[0], t[1], t[2]))

        f.write("#Part 3 - the hole list.\n")
        f.write('%d\n'%0)
        f.write("#Part 4 - the region list.\n")
        f.write('%d\n'%0)
'''

def surfgen_shared_boundary1(bbox=[-10,-10,-10, 10,10,10], l=3):
    # Generating point grids for two surfaces
    xmin = bbox[0]
    xmax  = bbox[3]
    ymin = bbox[1]
    ymax = bbox[4]
    zmin = bbox[2]
    zmax = bbox[5]
    X, Y, Z1 = sample_surf([xmax*0.8, ymax*0.8, zmax*0.3], step=0.8)
    Z2 = -Z1 - zmax*0.3
    Z1 = Z1 + zmax*0.3
    # Plotting the two surfaces
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(X, Y, Z1.reshape(-1,1), color='b')
    surf = ax.scatter(X, Y, Z2.reshape(-1,1), color='r')
    plt.show()
    
    # Integrating two surfaces
    points1 = np.concatenate((X.reshape(-1,1), Y.reshape(-1,1 ), Z1.reshape(-1,1)), axis=1)
    points2 = np.concatenate((X.reshape(-1,1), Y.reshape(-1,1 ), Z2.reshape(-1,1)), axis=1)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(points1[:,0],points1[:,1], points1[:,2], color='b')
    plt.show()
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(points2[:,0],points2[:,1], points2[:,2], color='b')
    plt.show()
    corners = utils.boundary_points(bbox)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(corners[:,0],corners[:,1], corners[:,2], color='b')
    plt.show()
    midcorners = utils.mid_corners(bbox)
    offset1 = len(corners) +len(midcorners) + 1
    offset2 = len(corners) + len(midcorners) + len(points1) + 1
    points = np.concatenate((corners, midcorners, points1, points2), axis=0)
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(points[:,0],points[:,1], points[:,2], color='b')
    plt.show()

    polys1 = list()
    polys2 = list()
    tmp1 = [9, 12]
    tmp2 = [9, 12]
    for x in np.sort(np.unique(X)):
        x_s1 = np.argwhere(points1[:,0]==x)
        x_s1_points = np.hstack((points1[x_s1.reshape(-1,)], x_s1))
        x_s1 = (x_s1_points[:,3][np.argsort(x_s1_points[:,1])] + offset1).astype(int)
        polys1.append(np.concatenate((tmp1[::-1], x_s1)))
        tmp1 = x_s1

        x_s2 = np.argwhere(points2[:,0]==x)
        x_s2_points = np.hstack((points2[x_s2.reshape(-1,)], x_s2))
        x_s2 = (x_s2_points[:,3][np.argsort(x_s2_points[:,1])] + offset2).astype(int)
        polys2.append(np.concatenate((tmp2[::-1], x_s2)))
        tmp2 = x_s2
    
    ymin_s1 = np.argwhere(points1[:,1] == np.amin(np.unique(Y)))
    ymin_s1_points = np.hstack((points1[ymin_s1.reshape(-1,)], ymin_s1))
    ymin_s1 = (ymin_s1_points[:,3][np.argsort(ymin_s1_points[:,0])] + offset1).astype(int)
    polys1.append(np.concatenate(([10, 9], ymin_s1)))

    ymax_s1 = np.argwhere(points1[:,1] == np.amax(np.unique(Y)))
    ymax_s1_points = np.hstack((points1[ymax_s1.reshape(-1,)], ymax_s1))
    ymax_s1 = (ymax_s1_points[:,3][np.argsort(ymax_s1_points[:,0])] + offset1).astype(int)
    polys1.append(np.concatenate(([11, 12], ymax_s1)))

    ymin_s2 = np.argwhere(points2[:,1] == np.amin(np.unique(Y)))
    ymin_s2_points = np.hstack((points2[ymin_s2.reshape(-1,)], ymin_s2))
    ymin_s2 = (ymin_s2_points[:,3][np.argsort(ymin_s2_points[:,0])] + offset2).astype(int)
    polys2.append(np.concatenate(([10, 9], ymin_s2)))

    ymax_s2 = np.argwhere(points2[:,1] == np.amax(np.unique(Y)))
    ymax_s2_points = np.hstack((points2[ymax_s2.reshape(-1,)], ymax_s2))
    ymax_s2 = (ymax_s2_points[:,3][np.argsort(ymax_s2_points[:,0])] + offset2).astype(int)
    polys2.append(np.concatenate(([11, 12], ymax_s2)))

    # Preparing PLC and save it to .poly file for tetgen
    with open('/home/altaa/tet/mediansurf.poly', 'w') as f:
        f.write("#Part 1 - the node list\n")
        f.write("#%d nodes in 3d, no attributes, no boundary marker\n"%points.shape[0])
        f.write('%d %d %d %d\n'%(points.shape[0], 3, 0,0))
        for i, p in enumerate(points):
            f.write("%d %f %f %f\n"%(i+1, p[0], p[1], p[2]))
        # Each 4 sides has 3 polygons
        # Top and bottom
        # Each triangle of the two surfaces are facets
        fn = 6 + len(polys1) + len(polys2)
        f.write("#Part 2 - the facet list.\n")
        f.write("#%d facets with boundary markers\n"%fn)
        f.write('%d %d\n'%(fn, 1))
        f.write("#Boundary facet list.\n")
        f.write("%d %d %d\n"%(1, 0, 1))
        f.write("4 1 2 3 4\n")
        f.write("%d %d %d\n"%(1, 0, 1))
        f.write("4 5 6 7 8\n")
        #xmin side
        f.write("2 0 1\n")
        f.write("4 1 4 8 5\n")
        f.write("2 9 12\n")
        #ymin side
        f.write("2 0 1\n")
        f.write("4 1 2 6 5\n")
        f.write("2 9 10\n")
        #xmax side
        f.write("2 0 1\n")
        f.write("4 2 3 7 6\n")
        f.write("2 10 11\n")
        #ymax side
        f.write("2 0 1\n")
        f.write("4 2 4 8 7\n")
        f.write("2 11 12\n")
        for p in polys1:
            f.write("%d %d %d\n"%(1, 0, -1))
            f.write("%d"%len(p) + " " + np.array_str(p)[1:-1] + "\n")
        for p in polys2:
            f.write("%d %d %d\n"%(1, 0, -2))
            f.write("%d"%len(p) + " " + np.array_str(p)[1:-1] + "\n")
        f.write("#Part 3 - the hole list.\n")
        f.write('%d\n'%0)
        f.write("#Part 4 - the region list.\n")
        f.write('%d\n'%0)

def surfgen(bbox=[-10,-10,-10, 10,10,10], l=2):
    # Generating point grids for two surfaces
    xmin = bbox[0]
    xmax  = bbox[3]
    ymin = bbox[1]
    ymax = bbox[4]
    zmin = bbox[2]
    zmax = bbox[5]
    mesh = meshgen2d([xmin, ymin, xmax, ymax], l, include_corners=True)
    #plot2d.plotmesh2d(mesh)
    #plt.show()
    X, Y, Z1 = sample_surf([xmax, ymax, zmax*0.3])
    Z2 = -Z1 - zmax*0.4
    Z1 = Z1 + zmax*0.4
    sample_points = np.hstack((X.reshape(-1,1), Y.reshape(-1,1)))
    # Plotting the two surfaces
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    surf = ax.scatter(X, Y, Z1.reshape(-1,1), color='b')
    surf = ax.scatter(X, Y, Z2.reshape(-1,1), color='r')
    plt.show()
    # Interpolating the surface mesh into two different surfaces
    # similar to the the sample surfaces generated before
    Z1 = interpolate_surf(sample_points, Z1.reshape(-1,1), mesh.points)
    Z2 = interpolate_surf(sample_points, Z2.reshape(-1,1), mesh.points)
    
    # Integrating two surfaces
    points1 = np.hstack((mesh.points, Z1))
    points2 = np.hstack((mesh.points, Z2))
    corners = utils.boundary_points(bbox)
    offset1 = len(corners) + 1
    offset2 = len(corners) + len(points1) + 1
    points = np.concatenate((corners, points1, points2), axis=0)

    triangles1 = mesh.simplices + offset1
    triangles2 = mesh.simplices + offset2
    triangles = np.vstack((mesh.simplices + offset1, mesh.simplices + offset2)) 
    # Plotting the two surface mesh
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(points1[:,0], points1[:,1], points1[:,2], triangles=mesh.simplices, cmap=cm.autumn, linewidth=0.2)
    ax.plot_trisurf(points2[:,0], points2[:,1], points2[:,2], triangles=mesh.simplices, cmap=cm.winter, linewidth=0.2)
    plt.show()

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_trisurf(points[:,0], points[:,1], points[:,2], triangles=triangles, cmap=cm.autumn, linewidth=0.2)
    surf = ax.scatter(points[:,0], points[:,1], points[:,2], color='y')
    plt.show()
    # 2 surfaces intersect with the 3D box and it results 8 curves on the 4 sides of a 3D box
    xmin_s1 = np.argwhere(points1[:,0]==xmin)
    xmin_s1_points = np.hstack((points1[xmin_s1.reshape(-1,)], xmin_s1))
    xmin_s1 = (xmin_s1_points[:,3][np.argsort(xmin_s1_points[:,1])] + offset1).astype(int)

    xmin_s2 = np.argwhere(points2[:,0]==xmin)
    xmin_s2_points = np.hstack((points2[xmin_s2.reshape(-1,)], xmin_s2))
    xmin_s2 = (xmin_s2_points[:,3][np.argsort(xmin_s2_points[:,1])] + offset2).astype(int)

    xmax_s1 = np.argwhere(points1[:,0]==xmax)
    xmax_s1_points = np.hstack((points1[xmax_s1.reshape(-1,)], xmax_s1))
    xmax_s1 = (xmax_s1_points[:,3][np.argsort(xmax_s1_points[:,1])] + offset1).astype(int)

    xmax_s2 = np.argwhere(points2[:,0]==xmax)
    xmax_s2_points = np.hstack((points2[xmax_s2.reshape(-1,)], xmax_s2))
    xmax_s2 = (xmax_s2_points[:,3][np.argsort(xmax_s2_points[:,1])] + offset2).astype(int)

    ymin_s1 = np.argwhere(points1[:,1]==ymin)
    ymin_s1_points = np.hstack((points1[ymin_s1.reshape(-1,)], ymin_s1))
    ymin_s1 = (ymin_s1_points[:,3][np.argsort(ymin_s1_points[:,0])] + offset1).astype(int)

    ymin_s2 = np.argwhere(points2[:,1]==ymin)
    ymin_s2_points = np.hstack((points2[ymin_s2.reshape(-1,)], ymin_s2))
    ymin_s2 = (ymin_s2_points[:,3][np.argsort(ymin_s2_points[:,0])] + offset2).astype(int)

    ymax_s1 = np.argwhere(points1[:,1]==ymax)
    ymax_s1_points = np.hstack((points1[ymax_s1.reshape(-1,)], ymax_s1))
    ymax_s1 = (ymax_s1_points[:,3][np.argsort(ymax_s1_points[:,0])] + offset1).astype(int)

    ymax_s2 = np.argwhere(points2[:,1]==ymax)
    ymax_s2_points = np.hstack((points2[ymax_s2.reshape(-1,)], ymax_s2))
    ymax_s2 = (ymax_s2_points[:,3][np.argsort(ymax_s2_points[:,0])] + offset2).astype(int)

    # Similarly each side of the box consists of 3 polygons
    polys = list()
    #xmin polygons
    polys.append(np.concatenate(([8, 5], xmin_s1)))
    polys.append(np.concatenate((xmin_s1[::-1], xmin_s2)))
    polys.append(np.concatenate((xmin_s2, [4, 1])))
    #xmax polygons
    polys.append(np.concatenate(([7, 6], xmax_s1)))
    polys.append(np.concatenate((xmax_s1[::-1], xmax_s2)))
    polys.append(np.concatenate((xmax_s2, [3, 2])))
    #ymin polygons
    polys.append(np.concatenate(([6, 5], ymin_s1)))
    polys.append(np.concatenate((ymin_s1[::-1], ymin_s2)))
    polys.append(np.concatenate((ymin_s2, [2, 1])))
    #ymax polygons
    polys.append(np.concatenate(([7, 8], ymax_s1)))
    polys.append(np.concatenate((ymax_s1[::-1], ymax_s2)))
    polys.append(np.concatenate((ymax_s2, [3, 4])))
    #zmin and zmax polygons
    polys.append(np.array([1, 2, 3, 4]))
    polys.append(np.array([5, 6, 7, 8]))
    # Preparing PLC and save it to .poly file for tetgen
    with open('/home/altaa/tet/tet.poly', 'w') as f:
        f.write("#Part 1 - the node list\n")
        f.write("#%d nodes in 3d, no attributes, no boundary marker\n"%points.shape[0])
        f.write('%d %d %d %d\n'%(points.shape[0], 3, 0,0))
        for i, p in enumerate(points):
            f.write("%d %f %f %f\n"%(i+1, p[0], p[1], p[2]))
        # Each 4 sides has 3 polygons
        # Top and bottom
        # Each triangle of the two surfaces are facets
        fn = 3*4 + 2 + len(triangles)
        f.write("#Part 2 - the facet list.\n")
        f.write("#%d facets with boundary markers\n"%fn)
        f.write('%d %d\n'%(fn, 1))
        f.write("#Boundary facet list.\n")
        for p in polys:
            f.write("%d %d %d\n"%(1, 0, 1))
            f.write("%d"%len(p) + " " + np.array_str(p)[1:-1] + "\n")
        f.write("#Facet list of surface1.\n")
        for t in triangles1:
            f.write("%d %d %d\n"%(1, 0, -1))
            f.write("%d %d %d %d\n"%(3, t[0], t[1], t[2]))
        f.write("#Facet list of surface2.\n")
        for t in triangles2:
            f.write("%d %d %d\n"%(1, 0, -2))
            f.write("%d %d %d %d\n"%(3, t[0], t[1], t[2]))
        f.write("#Part 3 - the hole list.\n")
        f.write('%d\n'%0)
        f.write("#Part 4 - the region list.\n")
        f.write('%d\n'%0)

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
'''
if __name__ == "__main__":
    surfgen_shared_boundary()
