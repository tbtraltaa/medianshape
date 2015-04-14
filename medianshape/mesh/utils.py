# encoding: utf-8

from __future__ import absolute_import

import sys
import numpy as np

from scipy.linalg import det
from scipy.misc import factorial
from scipy.spatial.distance import cdist, pdist
from scipy.sparse import dok_matrix

def extract_edges(simplices):
    edges = set()
    for simplex in simplices:
        for i in range(len(simplex)):
            edges.add(tuple(sorted([simplex[i], simplex[(i+1)%len(simplex)]])))
    return np.array(list(edges))

def get_subsimplices(simplices):
    simplices = np.sort(simplices, axis=1)
    subsimplices = set()
    for j in np.arange(simplices.shape[1]):
        idx = list(range(simplices.shape[1]))
        idx.pop(j)
        subsimplices = subsimplices.union(set(tuple(sub) for sub in simplices.take(idx, axis=1)))
    return np.array([ sub for sub in subsimplices])

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
    return simplices

def simpvol(points, simplices):

    ''' SIMPVOL Simplex volume.
        V=SIMPVOL(points, simplices)
        Copyright (C) 2004-2006 Per-Olof Persson. See COPYRIGHT.TXT for details.
    '''
    point_dim  =  points.shape[1]
    simp_dim = simplices.shape[1]
    # 1-simplex, edge in 1 dimension
    if simp_dim == 1:
        volume = np.ones(simplices.shape[0], dtype=int) 
    elif point_dim  ==  1 and simp_dim == 2:
        d12 = points[simplices[:,1],:] - points[simplices[:,0],:]
        volume = np.abs(d12)
    # 1-simplex, edge in any dimension
    elif simp_dim == 2:
        d12 = points[simplices[:,1],:] - points[simplices[:,0],:]
        volume = np.sqrt(np.sum(d12**2, axis=1))
    # 2-simplex, triangle in 2 dimension
    elif point_dim  ==  2 and simp_dim == 3:
        d12 = points[simplices[:,1],:] - points[simplices[:,0],:]
        d13 = points[simplices[:,2],:] - points[simplices[:,0],:]
        volume = (np.multiply(d12[:,0], d13[:,1]) - np.multiply(d12[:,1], d13[:,0]))*1.0/2
    # 2-simplex, triangle in 3 dimension
    elif point_dim == 3 and simp_dim == 3:
        d12 = points[simplices[:,1],:] - points[simplices[:,0],:]
        d13 = points[simplices[:,2],:] - points[simplices[:,0],:]
        volume = np.sqrt(np.sum(np.cross(d12, d13, axis=1)**2, axis=1))*1.0/2
    # 3-simplex, tetrahedra in 3 dimension
    elif point_dim  ==  3 and simp_dim == 4: # 
        d12 = points[simplices[:,1],:] - points[simplices[:,0],:]
        d13 = points[simplices[:,2],:] - points[simplices[:,0],:]
        d14 = points[simplices[:,3],:] - points[simplices[:,0],:]
        volume = np.dot(np.cross(d12,d13,axis=2),d14,axis=2)*1.0/6
    # n-simplex in n-dimention
    else:
        volume = np.zeros((simplices.shape[0],1))
        for i in range(simplices.shape[0]):
            A = np.zeros(points.shape[1] + 1)
            A[:,0] = 1
            for j in range(points.shape[1]+1):
                A[j,1:] = points[simplices[i,j],:]
            volume[i] = np.det(A)
        volume = volume/np.factorial(points.shape[1])
    return volume

def simpvol1(points, simplices):
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
    if dim == 2 and n == 1:
        volume = 1
    elif dim == 2 and n == 3:
        for i, simplex in enumerate(simplices):
            extended_simplex = points[simplex]
            extended_simplex = np.hstack((extended_simplex, np.ones((n,1))))
            volume[i] = det(extended_simplex)/factorial(dim)
    elif dim == 2 and n == 2:
        for i, simplex in enumerate(simplices):
            volume[i] = pdist(points[simplex])
    return volume
    

# Builds a boundary matrix of given simplices. The format of a boundary matrix is as follows.
# boundary_matrix = (number of edges) x (number of simplices)
def boundary_matrix(simplices, subsimplices, is_oriented=True, is_sparse=True, format='coo'):
    simplex_dim  = simplices.shape[1] 
    #if simplex_dim - subsimplices.shape[1] != 1:
    #    sys.stderr.write("Unable to build a boundary matrix. Please enter (d+1)-simplices and  d-subsimplices\n")
     #   exit()
    n_simplices = simplices.shape[0]
    m_subsimplices = subsimplices.shape[0] 
    if is_sparse:
        boundary_matrix = dok_matrix((m_subsimplices, n_simplices), dtype=np.int8) 
    else:
        boundary_matrix = np.array((m_subsimplices, n_simplices), dtype=np.int8) 
    if is_oriented:
        simplices_sort_idx = np.argsort(simplices)
        subsimplices_sort_idx = np.argsort(subsimplices)
        simplices_parity = permutationparity(simplices_sort_idx, 2)
        subsimplices_parity = permutationparity(subsimplices_sort_idx, 2)
    val = 1
    simplices = np.sort(simplices, axis=1)
    subsimplices = np.sort(subsimplices, axis=1)
    for i, simplex in enumerate(simplices):
        for j in np.arange(simplex_dim):
            idx = list(range(simplex_dim))
            idx.pop(j)
            subsimplex = simplex.take(idx)
            # to check the membership of subsimplex in subsimplices
            subsimplex_idx = np.argwhere((subsimplices==subsimplex).all(axis=1) == True)
            if subsimplex_idx.size == 0:
                sys.stderr.write("Unable to find subsimplex! Make sure subsimplices contains \
                all boundary subsimplices")
                exit()
            subsimplex_idx = subsimplex_idx[0][0]
            if is_oriented:
                val = (-1)**((j + 1) + 1 + simplices_parity[i] + subsimplices_parity[subsimplex_idx])
            boundary_matrix[subsimplex_idx, i] = val
    if is_sparse:
        return boundary_matrix.asformat(format)
    else:
        return boundary_matrix

def boundary_matrix_2d(simplices, edges, sparse=True, format=None):
    if not sparse:
        boundary_matrix = np.array((len(edges), len(simplices)), dtype=np.int8) 
    else:
        boundary_matrix = dok_matrix((len(edges), len(simplices)), dtype=np.int8) 
    temp_edges = [ tuple(edge) for edge in edges]
    for j, simplex in enumerate(simplices):
        if len(simplex) == 2:
            boundary_matrix[simplex[0],j] = 1
            boundary_matrix[simplex[1],j] = -1
        else:
            simplex_boundary = boundary(simplex)
            for s_edge in simplex_boundary:
                sorted_s_edge = np.sort(s_edge)
                i = temp_edges.index(tuple(sorted_s_edge)) 
                if s_edge[0] <= s_edge[1]:
                    boundary_matrix[i,j] = 1
                else:
                    boundary_matrix[i,j] = -1
    if sparse and format:
        return boundary_matrix.asformat(format)
    else:
        return boundary_matrix

# Returns simplex boundary as faces.
# if an index given, returns n-1 simplex by removing the element at the index.
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

def permutationparity(P, dim=None):
    nRows, mCols = P.shape
    if nRows == 1 and dim is None:
        dim = 2
    elif mCols == 1 and dim is None:
        dim = 1
    p = [0]
    if dim == 1:
        for i in np.arange(nRows):
            p = np.sum(np.tile(P[i,:].reshape(-1,1), (nRows-(i+1), 1)) > P[i+1:,:], 0) + p
    elif dim == 2:
        for j in np.arange(mCols):
            p = np.sum(np.tile(P[:, j].reshape(-1,1), (1, mCols-(j+1))) > P[:, j+1:], 1) + p
    p = np.mod(p, 2)
    return p

if __name__ == '__main__':
    points = np.array([[0, 0, 0],[0,1,1],[0,2,0], [2,2,0]])
    simplices = np.array([[0, 1, 2], [1,2,3], [2, 4, 3]])
    simplices = np.array([[2, 1, 0], [3,1,2], [3, 2, 4]])
    permutationparity(simplices)
    simpvol1(points, simplices)

