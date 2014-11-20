# encoding: utf-8

from __future__ import absolute_import

import numpy as np

from scipy.linalg import det
from scipy.misc import factorial
from scipy.spatial.distance import cdist, pdist

def extract_edges(simplices):
    edges = set()
    for simplex in simplices:
        for i in range(len(simplex)):
            edges.add(tuple(sorted([simplex[i], simplex[(i+1)%len(simplex)]])))
    return np.array(list(edges))

def orient_simplices(simplices):
    direction = right_hand_rule(simplices[0])
    simplices = np.array(range(0, len(simplices)))
    if direction < 0:
        simplices[0] = simplices[0,::-1]
    for i, simplex in enumerate(simplices):
        simplices = np.delete(simplices, np.where(simplices==i))
        neighbors = self.mesh.neighbors[i]
        for opposit_point in np.where(neighbors >= 0)[0]:
            n_simplex_idx = neighbors[opposit_point]
            if any(simplices==n_simplex_idx):
                n_simplex = simplices[n_simplex_idx]
                n_boundary = Mesh.boundary(n_simplex)
                subsimplex= Mesh.boundary(simplex, opposit_point)
                for n_face in n_boundary:
                    if all((np.array(subsimplex) - np.array(n_face)) == 0):
                        simplices[n_simplex_idx] = n_simplex[::-1]
                simplices = np.delete(simplices, np.where(simplices==n_simplex_idx))
    print right_hand_rule(simplices[0])
    return simplices

def simpvol(points, simplices):
    n = simplices.shape[1]
    dim = points.shape[1]
    volume = np.zeros(simplices.shape[0])
#    if dim == 1:
#        d01 = p[t[:,1]]-p[t[:,0]]
#        return d01
#    elif dim == 2:
#        d01 = p[t[:,1]]-p[t[:,0]]
#        d02 = p[t[:,2]]-p[t[:,0]]
#        return (d01[:,0]*d02[:,1]-d01[:,1]*d02[:,0])/2
    if dim == 2 and n == 3:
        for i, simplex in enumerate(simplices):
            extended_simplex = points[simplex]
            extended_simplex = np.hstack((extended_simplex, np.ones((n,1))))
            volume[i] = det(extended_simplex)/factorial(n)
    elif dim == 2 and n == 2:
        for i, simplex in enumerate(simplices):
            volume[i] = pdist(points[simplex])
    return volume
    
def boundary1(simplex):
    boundary = list()
    if len(simplex) > 2:
        n = len(simplex)
        boundary.append(simplex[1:])
        for i in xrange(1,n):
            face = np.append(simplex[:i], simplex[i+1:])
            if (-1)^(i+1) < 0: 
                face = face[::-1]
            boundary.append(face)
    return boundary

def boundary_matrix(simplices, edges):
    boundary_matrix = np.zeros((len(edges), len(simplices))) 
    for i, edge in enumerate(edges):
        for j, simplex in enumerate(simplices):
            simplex_boundary = boundary(simplex)
            for s_edge in simplex_boundary:
                if all((s_edge - edge) == 0):
                    boundary_matrix[i,j] = 1
                    break
                elif all((s_edge - np.array(list(reversed(edge)))) == 0):
                    boundary_matrix[i,j] = -1
                    break
    return boundary_matrix

def boundary(simplex, idx=None):
    boundary = list()
    if idx == None:
        n = len(simplex)
        for i in xrange(0,n):
            face = list(simplex)
            face.pop(i)
            if (-1)**(i+2) < 0: 
                face = face[::-1]
            boundary.append(face)
    else:
        face = list(simplex)
        face.pop(idx)
        if (-1)**(idx+2) < 0: 
            face = face[::-1]
        boundary = face
    return boundary