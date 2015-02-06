# encoding: utf-8

from __future__ import division

import importlib
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, cdist

import utils

def generate_curve_on_mesh(points, mesh, interval_size=30, is_closed=False, func_str=None):
        closest_vertices = find_closest_vertices(points, mesh, interval_size)
        print "closest_vertices:\n", closest_vertices
        print "Function points:\n", points
        curve_path = find_path(closest_vertices, mesh, is_closed)
        #Mesh.disp_edges(func_path, "Function edges")
        edge_vector = get_edge_vector(curve_path, mesh)
        return edge_vector, curve_path, closest_vertices

def find_closest_vertices(points, mesh, interval_size=1, func_str=None):
    closest_vertices = list()
    for point in points:
        closest_vertex = find_closest_vertex(point, mesh, closest_vertices, interval_size, func_str)
        closest_vertices.append(closest_vertex)
    return  np.array(closest_vertices)

def find_closest_vertex(point, mesh, selected_points, interval_size=5, func_str=None):
    ordered_X = np.sort(np.unique(mesh.points[:,0]))
    interval_X = find_interval_X(point, ordered_X, interval_size)
    min_dist = 1000
    closest_vertex = -1
    min_idx = 0
    if func_str:
        func_values = utils.vectorize(func_str, interval_X)
        func_points = np.hstack((interval_X, func_values))
    for i, x in enumerate(interval_X):
        #Assigning candidate indices out of a tuple with 2 elements where the second one is empty
        candidate_idx = np.where(mesh.points[:,0]==x)[0]     
        candidate_idx = [c_idx for i, c_idx in enumerate(candidate_idx) if c_idx not in selected_points]
        if len(candidate_idx) != 0:
            candidate_points = mesh.points[candidate_idx]
            if func_str:
                distances = cdist(candidate_points, func_points[i].T.reshape(1,2))
            else:
                distances = cdist(candidate_points, point.T.reshape(1,2))
            min_idx = np.argmin(distances)  
            dist = distances[min_idx]
            if min_idx.size > 1:
                min_idx = random.sample(min_idx, 1)
            if dist < min_dist:
                min_dist = dist
                closest_vertex = candidate_idx[min_idx]
    return closest_vertex


def find_interval_X(point, ordered_X, interval_size=5):
    x_idx = np.where(ordered_X==point[0])
    if x_idx [0]:
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

def find_path(path_points, mesh, is_closed=False):
    no_of_points = mesh.points.shape[0] 
    adjacency_matrix = np.zeros((no_of_points, no_of_points), dtype=np.float64)  
    if is_closed:
        path_points = np.append(path_points, path_points[0].reshape(1,), axis=0)
    for i, edge in enumerate(mesh.edges):
            p1 = mesh.points[edge[0]].T.reshape(1,2)
            p2 = mesh.points[edge[1]].T.reshape(1,2)
            adjacency_matrix[edge[0],edge[1]] = cdist(p1, p2, 'euclidean')[0]

    #graph = csr_matrix(adjacency_matrix)  
    graph = adjacency_matrix
    #print "Adjacent matrix \n", adjacency_matrix
    #print graph.data

    # Finding path vertices in reverse order from the end to the start
    path_vertices = list()
    for i, point in reversed(list(enumerate(path_points))): 
        if i-1 >= 0:
            i1 = path_points[i-1]
            i2 = point
            distances, predecessors = dijkstra(graph, directed=False, \
                                        indices=i1, return_predecessors=True)
            j = i2
            while j != i1:
                if len(path_vertices) == 0:
                    path_vertices.append(j)
                    graph[i][j] = 100
                    #graph[j][1] = graph[j][0][1] + 100
                    #print graph[j]

                elif len(path_vertices) > 0 and path_vertices[-1]!= j:
                    path_vertices.append(j)
                    graph[i][j] = 100
                    #print graph[j]
                    #graph[j][1] = graph[j][0][1] + 100
                j = predecessors[j]
            path_vertices.append(i1)
    path = list()
    print "Path", path_vertices
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

def get_edge_vector(path, mesh):
    edge_vector = np.zeros(shape=(mesh.edges.shape[0], 1))
    for i, edge in enumerate(mesh.edges):
        for path_edge in path:
            if all((path_edge - edge) == 0):
                edge_vector[i] = 1
            elif all((path_edge - np.array(list(reversed(edge)))) == 0):
                edge_vector[i] = -1
    return edge_vector

def plot_curve(mesh, input_points, closest_vertices, path, title=None, color="red", linewidth=3, label=""):
    plt.plot(input_points[:,0], input_points[:,1], c=color, ls="--")
    plt.title(title)
    for i, edge in enumerate(path):
        points = mesh.points[edge]
        plt.plot(points[:,0], points[:,1], color, linewidth=linewidth, label=label)
    plt.scatter(mesh.points[closest_vertices][:,0], mesh.points[closest_vertices][:,1], s=100)
    plt.scatter(input_points[:,0], input_points[:,1], c=color)

if __name__ == '__main__':
    pass
