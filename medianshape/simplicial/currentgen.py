# encoding: utf-8
'''
Current generation
==================

'''

from __future__ import division

import sys
import os
import importlib

import numpy as np
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import dok_matrix
from scipy.spatial.distance import pdist, cdist

from medianshape.simplicial.utils import vectorize

def push_functions_on_mesh_2d(mesh, curves, functions, is_closed=False):
    '''
    Hi
    '''
    input_currents = list()
    paths = list()
    vertices = list()
    for i, curve_points in enumerate(curves):
        func_str = functions[i]
        input_current, path, closest_vertices = \
        push_function_on_mesh(mesh, curve_points, func_str=func_str, is_closed=is_closed) 
        vertices.append(closest_vertices)
        paths.append(path)
        input_currents.append(input_current)
    input_currents = np.array(input_currents).reshape(len(curves), mesh.edges.shape[0])
    return vertices, paths, input_currents

def push_function_on_mesh(mesh, points, interval_size=10, func_str=None, is_closed=False):
    '''
    HI
    '''
    closest_vertices = list()
    for point in points:
        min_dist = mesh.diagonal
        closest_vertex = -1
        min_idx = 0
        ordered_X = np.sort(np.unique(mesh.points[:,0]))
        interval_X = find_interval_X(point, ordered_X, interval_size)
        func_values = vectorize(func_str, interval_X)
        func_points = np.hstack((interval_X, func_values))
        for i, x in enumerate(interval_X):
            #Assigning candidate indices out of a tuple with 2 elements where the second one is empty
            candidate_idx = np.where(mesh.points[:,0]==x)[0]     
            candidate_idx = [c_idx for c_idx in candidate_idx if c_idx not in closest_vertices]
            if len(candidate_idx) != 0:
                candidate_points = mesh.points[candidate_idx]
                distances = cdist(candidate_points, func_points[i].reshape(1,-1))
                min_idx = np.argmin(distances)  
                dist = distances[min_idx]
                if min_idx.size > 1:
                    min_idx = random.sample(min_idx, 1)
                if dist < min_dist:
                    min_dist = dist
                    closest_vertex = candidate_idx[min_idx]
        closest_vertices.append(closest_vertex)
    closest_vertices = np.array(closest_vertices)
    curve_path = find_path(mesh, closest_vertices, is_closed)
    edge_vector = get_edge_vector(mesh, curve_path)
    return edge_vector, curve_path, closest_vertices

def push_curves_on_mesh(mesh, simplices, subsimplices, curves, is_closed=False, valid_points=None):
    '''
    Hi
    '''
    input_currents = list()
    paths = list()
    vertices = list()
    for i, curve_points in enumerate(curves):
        input_current, path, closest_vertices = \
        push_curve_on_mesh(mesh, curve_points, is_closed, valid_points) 
        vertices.append(closest_vertices)
        paths.append(path)
        input_currents.append(input_current)
    input_currents = np.array(input_currents).reshape(len(curves), subsimplices.shape[0])
    return vertices, paths, input_currents

def push_curve_on_mesh(mesh, points, is_closed=False, valid_points=None):
    '''
    Hi
    '''
    closest_vertices = find_closest_vertices(mesh, points, valid_points)
    if len(closest_vertices) == 1:
        input_current = get_vertex_vector(mesh, closest_vertices)
        curve_path = closest_vertices
    else:
        curve_path = find_path(mesh, closest_vertices, is_closed)
        input_current = get_edge_vector(mesh, curve_path)
    return input_current, curve_path, closest_vertices

def find_closest_vertices(mesh, points, valid_points=None):
    '''
    Hi
    '''
    closest_vertices = list()
    if valid_points is not None and len(points) > len(valid_points):
        sys.stderr.write("Too many points given.\n")
        exit()
    if valid_points is None and len(points) > len(mesh.points):
        sys.stderr.write("Number of input points are more than the number of points in the mesh.\n")
        exit()
    for point in points:
        closest_vertex = find_closest_vertex(mesh, point, closest_vertices, valid_points)
        closest_vertices.append(closest_vertex)
    return  np.array(closest_vertices)

def find_closest_vertex(mesh, point, selected_points, valid_points=None):
    '''
    Hi
    '''
    min_dist = mesh.diagonal
    closest_vertex = -1
    min_idx = 0
    temp_points = np.copy(mesh.points)
    distances = cdist(point.reshape(1, -1), temp_points)
    closest_vertex = np.argmin(distances)  
    while closest_vertex in selected_points:
        temp_points[closest_vertex] = temp_points[np.argmax(distances)]
        distances = cdist(point.reshape(1, -1), temp_points)
        closest_vertex= np.argmin(distances)  
    if valid_points is not None:
        while closest_vertex not in valid_points or closest_vertex in selected_points:
            temp_points[closest_vertex] = temp_points[np.argmax(distances)]
            distances = cdist(point.reshape(1, -1), temp_points)
            closest_vertex= np.argmin(distances)  
    return closest_vertex

def find_interval_X(point, ordered_X, interval_size=5):
    '''
    Hi
    '''
    x_idx = np.where(ordered_X==point[0])
    if len(x_idx[0]):
        x_idx = ordered_X.tolist().index(point[0])
    else:
        next_x = [x for x in ordered_X if x >= point[0]][0]
        x_idx = ordered_X.tolist().index(next_x)
    interval_X = []      
    i = 0
    start_idx = 0
    while True:
        interval_idx = x_idx - interval_size//2 + i
        if interval_idx >= 0:
            start_idx = interval_idx
            break
        i += 1
    next_x_idx = x_idx + interval_size
    end_idx = start_idx + interval_size
    next_start_idx = next_x_idx - interval_size//2
    if end_idx >= next_start_idx:
        interval_X = ordered_X[start_idx:next_start_idx]
    else:
        interval_X = ordered_X[start_idx:start_idx + interval_size]
    return interval_X.reshape(interval_X.size, 1)

def find_path(mesh, path_points, is_closed=False):
    '''
    Hi
    '''
    no_of_points = mesh.points.shape[0] 
    adjacency_matrix = dok_matrix((no_of_points, no_of_points), dtype=np.float64)  
    if is_closed:
        path_points = np.append(path_points, path_points[0].reshape(1,), axis=0)
    for i, edge in enumerate(mesh.edges):
            p1 = mesh.points[edge[0]].T.reshape(1,-1)
            p2 = mesh.points[edge[1]].T.reshape(1,-1)
            adjacency_matrix[edge[0],edge[1]] = cdist(p1, p2, 'euclidean')[0]
            adjacency_matrix[edge[1],edge[0]] = cdist(p1, p2, 'euclidean')[0]

    # Finding path vertices in reverse order from the end to the start
    path_vertices = list()
    for i, point in reversed(list(enumerate(path_points))): 
        if i-1 >= 0:
            i1 = path_points[i-1]
            i2 = point
            distances, predecessors = dijkstra(adjacency_matrix, directed=False, \
                                        indices=i1, return_predecessors=True)
            j = i2
            while j != i1:
                if len(path_vertices) == 0:
                    path_vertices.append(j)

                elif len(path_vertices) > 0 and path_vertices[-1]!= j:
                    path_vertices.append(j)
                    adjacency_matrix[prev,j] = mesh.diagonal
                    adjacency_matrix[j,prev] = mesh.diagonal
                prev = j    
                j = predecessors[j]
            adjacency_matrix[prev,i1] = mesh.diagonal
            adjacency_matrix[i1,prev] = mesh.diagonal
            prev = i1
            path_vertices.append(i1)

    path = list()
    # Generating path made of edges based on the path edges
    # The path edges are not reverse order
    for i, point in enumerate(path_vertices):
        if i+1 < len(path_vertices):
            edge = list([path_vertices[i+1], point])
            path.insert(0, edge)
    return  np.array(path, dtype=int)

#    def find_neighbors(self, point_idx):
#        return self.mesh.vertex_neighbor_vertices[1][self.mesh.vertex_neighbor_vertices[0]\
#        [point_idx]:self.mesh.vertex_neighbor_vertices[0][point_idx+1]]

def get_edge_vector(mesh, path):
    '''
    Hi
    '''
    edge_vector = np.zeros(shape=(mesh.edges.shape[0], 1))
    temp_edges = [ tuple(edge) for edge in mesh.edges]
    for path_edge in path:
        sorted_path_edge = np.sort(path_edge)
        i = temp_edges.index(tuple(sorted_path_edge))
        if path_edge[0] <= path_edge[1]:
            edge_vector[i] = 1
        else:
            edge_vector[i] = -1
    return edge_vector

def get_vertex_vector(mesh, path):
    '''
    Hi
    '''
    vertex_vector = np.zeros(shape=(mesh.points.shape[0], 1))
    vertex_vector[path] = 1
    return vertex_vector
            
if __name__ == '__main__':
    pass
