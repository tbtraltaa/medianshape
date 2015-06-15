# encoding: utf-8

from __future__ import absolute_import

import importlib

import numpy as np

import distmesh as dm
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt

import medianshape.utils as utils 
from medianshape.simplicial.mesh import Mesh2D, Mesh3D

def meshgen2d(boundary_box=None, l=0.02, fixed_points=None, include_corners=True):
    mesh = Mesh2D()
    #l - initial length of triangle sides. Change it to vary traingle size
    mesh.bbox = boundary_box
    mesh.set_boundary_points()
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.fixed_points = fixed_points
    if include_corners:
        if mesh.fixed_points is not None:
            mesh.fixed_points = np.vstack((mesh.fixed_points, mesh.boundary_points))
        else:
            mesh.fixed_points = mesh.boundary_points
    mesh.points, mesh.simplices = distmesh2d('square', mesh.bbox, mesh.fixed_points, l=l)
    mesh.edges = utils.get_subsimplices(mesh.simplices)
    mesh.orient_simplices_2D()
    return mesh

def meshgen3d(boundary_box=None, l=0.2, fixed_points=None, include_corners=True, load_data=False):
    mesh = Mesh3D()
    mesh.bbox = boundary_box
    mesh.set_boundary_points()
    mesh.set_diagonal()
    mesh.set_boundary_values()
    mesh.fixed_points = fixed_points
    if include_corners:
        if mesh.fixed_points is not None:
            mesh.fixed_points = np.vstack((mesh.fixed_points, mesh.boundary_points))
        else:
            mesh.fixed_points = mesh.boundary_points
    mesh.points, mesh.simplices= distmesh3d("sphere", mesh.bbox, mesh.fixed_points, l=l)
    #mesh.points, mesh.simplices = scipy_mesh3d(mesh.bbox, mesh.fixed_points, l)
    mesh.triangles = utils.get_subsimplices(mesh.simplices)
    mesh.edges = utils.get_subsimplices(mesh.triangles)
    return mesh

def scipy_mesh3d(bbox, fixed_points, l):
    bbox = np.array(bbox).reshape(2, -1)
    dim = bbox.shape[1]
    points = np.mgrid[tuple(slice(min, max+l, l) for min, max in bbox.T)]
    points = points.reshape(dim, -1).T
    ax = plt.gca(projection='3d')
    ax.scatter(points[:,0],points[:,1], points[:,2])
    return points, Delaunay(points).simplices

#def meshpy_cuboid3d(bbox, fixed_points, max_volume):
#    mesh_info = MeshInfo()
#    mesh_info.set_points(fixed_points)
#    mesh_info.set_facets([[0,1,2,3],\
#                            [4,5,6,7],\
#                            [0,4,5,1],\
#                            [1,5,6,2],\
#                            [2,6,7,3],\
#                            [3,7,4,0],\
#                                ])
#    mesh = build(mesh_info, volume_constraints=True, max_volume=max_volume)
#    return np.array(mesh.points), np.array(mesh.elements)

def distmesh2d(shape, bbox, fixed_points, l=0.1):
    if shape == "square":
        return square_mesh(bbox, fixed_points, l)

def distmesh3d(shape, bbox, fixed_points, l=0.1):
    #if shape == "cuboid":
        #return cuboid_mesh(bbox, fixed_points, l)
    if shape == "sphere":
        return sphere_mesh(bbox, fixed_points, l=l)

def square_mesh(bbox, fixed_points, l=0.1):
    """Square, with size function point and line sources"""
    dist_function = lambda p: dm.drectangle(p, bbox[0], bbox[2], \
                                                bbox[1],bbox[3])
    points, triangles = dm.distmesh2d(dist_function, dm.huniform, l, bbox, fixed_points)
    return points, triangles 

def cuboid_mesh(bbox, fixed_points, l=0.1):
    '''
        Generates 3-simplex, tetrahedral mesh in 3D using DistMesh
    '''
    dist_function = lambda p: dcuboid(p, float(bbox[0]), float(bbox[1]), float(bbox[2]), float(bbox[3]), float(bbox[4]), float(bbox[5])) 
    points, tetrahedras = dm.distmeshnd(dist_function, dm.huniform, l, bbox, fixed_points) 
    return points, tetrahedras

def dcuboid(p, x1, y1, z1, x2, y2, z2):
    """Signed distance function for cubiod with corners (x1,y1,z1), (x1,y1,z2),
    (x1,y2,z1), (x1,y2,z2), (x2,y1,z1), (x2,y2,z1), (x2,y1,z2), (x2,y2,z2).

    This has an incorrect distance to the four corners. 
    """
    d1 = x1 - p[:,0]
    d2 = -x2 + p[:,0]
    d3 = y1 - p[:,1]
    d4 = -y2  + p[:,1]
    d5 = z1 - p[:,2]
    d6 = -z2 + p[:,2]

    d135 = np.sqrt(d1**2+d3**2+d5**2)
    d136 = np.sqrt(d1**2+d3**2+d6**2)
    d145 = np.sqrt(d1**2+d4**2+d5**2)
    d146 = np.sqrt(d1**2+d4**2+d6**2)
    d235 = np.sqrt(d2**2+d3**2+d5**2)
    d236 = np.sqrt(d2**2+d3**2+d6**2)
    d245 = np.sqrt(d2**2+d4**2+d5**2)
    d246 = np.sqrt(d2**2+d4**2+d6**2)
    
    d = -np.minimum(np.minimum(np.minimum(np.minimum(np.minimum(-d1,-d2),-d3),-d4),-d5),-d6)
    ix = (d1>0)*(d3>0)*(d5>0)
    d[ix] = d135[ix]
    ix = (d1>0)*(d3>0)*(d6>0)
    d[ix] = d136[ix]
    ix = (d1>0)*(d4>0)*(d5>0)
    d[ix] = d145[ix]
    ix = (d1>0)*(d4>0)*(d6>0)
    d[ix] = d146[ix]
    ix = (d2>0)*(d3>0)*(d5>0)
    d[ix] = d235[ix]
    ix = (d2>0)*(d3>0)*(d6>0)
    d[ix] = d236[ix]
    ix = (d2>0)*(d4>0)*(d5>0)
    d[ix] = d245[ix]
    ix = (d2>0)*(d4>0)*(d4>0)
    d[ix] = d246[ix]
    return d

def sphere_mesh(bbox, fixed_points, l):
    r = np.abs(bbox[3] - bbox[0])*1.0/2  
    center = [(bbox[0]+bbox[3])*1.0/2, (bbox[1] + bbox[4])*1.0/2, (bbox[2]+bbox[5])*1.0/2]
    dist_function = lambda p: dm.dsphere(p, center[0], center[1], center[2], r)
    points, simplices = dm.distmeshnd(dist_function, dm.huniform, l, bbox, fixed_points, fig=None) 
    return points, simplices

