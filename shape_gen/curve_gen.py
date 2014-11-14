# encoding: utf-8

from __future__ import division, absolute_import

import importlib
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix
from scipy.spatial.distance import pdist, cdist


class FunctionApprox2d():
    def __init__(self, mesh, func_str=None, sample_step=None, interval_size=None):
        self.mesh = mesh
        self.points = mesh.points
        self.simplices = mesh.simplices
        self.edges = mesh.edges
        self.func_str = func_str
        self.ordered_X = np.sort(np.unique(self.points[:,0]))

        sample_size = int(math.sqrt(len(self.points)))
        if not sample_step:
            self.sample_step = self.ordered_X.size//sample_size
        if interval_size and interval_size < sample_step -1:
            self.interval_size = interval_size
        else:
            self.interval_size = self.sample_step

    def generate_curve(self, func_str=None):
        if func_str:
            self.func_str = func_str
        self.nearest_points, self.func_points = self.find_nearest_points()
        print "Nearest_points:\n", self.nearest_points
        print "Function points:\n", self.func_points
        self.func_path = self.find_path(self.nearest_points)
        #Mesh.disp_edges(func_path, "Function edges")
        path_vector = self.get_path_vector()
        return path_vector

    def find_nearest_points(self):
        sample_X = []
        nearest_points = list()
        for i in range(0, self.ordered_X.size, self.sample_step):
            nearest_points.append(self.find_nearest_point(i))
            sample_X.append(self.ordered_X[i])
        func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), \
        FunctionApprox2d.vectorize(self.func_str, sample_X).reshape(len(sample_X),1)))
        return  np.array(nearest_points), func_points

    def find_nearest_point(self, x_idx):
        interval_X = self.find_interval_X(x_idx)
        func_values = FunctionApprox2d.vectorize(self.func_str, interval_X)
        min_dist = 1000
        nearest_point = -1
        func_points = np.hstack((interval_X, func_values))
        for i, x in enumerate(interval_X):
            #Assigning candidate indices out of a tuple with 2 elements where the second one is empty
            candidate_idx = np.where(self.points[:,0]==x)[0]     
            candidate_points = self.points[candidate_idx]
            distances = cdist(candidate_points, func_points[i].T.reshape(1,2))
            min_idx = np.argmin(distances)  
            dist = distances[min_idx]
            if min_idx.size > 1:
                min_idx = random.sample(min_idx, 1)
            if dist < min_dist:
                min_dist = dist
                nearest_point = candidate_idx[min_idx]
        return nearest_point

    def find_interval_X(self, x_idx):
        interval_X = []      
        i = 0
        start_idx = 0
        while True:
            interval_idx = x_idx - self.interval_size//2 + i
            if interval_idx >= 0:
                start_idx = interval_idx
                break
            i += 1
        next_x_idx = x_idx + self.interval_size
        end_idx = start_idx + self.interval_size
        next_start_idx = next_x_idx - self.interval_size//2
        if end_idx >= next_start_idx:
            interval_X = self.ordered_X[start_idx:next_start_idx]
        else:
            interval_X = self.ordered_X[start_idx:start_idx + self.interval_size]
        return interval_X.reshape(interval_X.size, 1)

    def find_path(self, nearest_points):
        no_of_points = self.points.shape[0] 
        adjacency_matrix = np.zeros((no_of_points, no_of_points), dtype=np.float64)  
        for i, edge in enumerate(self.edges):
                p1 = self.points[edge[0]].T.reshape(1,2)
                p2 = self.points[edge[1]].T.reshape(1,2)
                adjacency_matrix[edge[0],edge[1]] = cdist(p1, p2, 'euclidean')[0]

        #graph = csr_matrix(adjacency_matrix)  
        graph = adjacency_matrix
        #print "Adjacent matrix \n", adjacency_matrix
        #print graph.data
        path_vertices = list()
        for i, point in reversed(list(enumerate(nearest_points))): 
            if i-1 >= 0:
                i1 = nearest_points[i-1]
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
        path_vertices = np.array(path_vertices)
        print "Path", path_vertices
        for i, point in reversed(list(enumerate(path_vertices))):
            if i+1 < len(path_vertices):
                edge = list([path_vertices[i+1], point])
                path.append(edge)
        return  np.array(path, dtype=int)

#    def find_neighbors(self, point_idx):
#        return self.mesh.vertex_neighbor_vertices[1][self.mesh.vertex_neighbor_vertices[0]\
#        [point_idx]:self.mesh.vertex_neighbor_vertices[0][point_idx+1]]

    def get_path_vector(self):
        path_vector = np.zeros(shape=(self.edges.shape[0], 1))
        for i, edge in enumerate(self.edges):
            for path_edge in self.func_path:
                if all((path_edge - edge) == 0):
                    path_vector[i] = 1
                    break
                elif all((path_edge - np.array(list(reversed(edge)))) == 0):
                    path_vector[i] = -1
                    break
        return path_vector
    
    @staticmethod
    def vectorize(func_str, X):
        func_points = []
        if func_str.find(".") != -1:
            mod_name, func_name = func_str.rsplit('.', 1)
            mod = importlib.import_module(mod_name)
            func = getattr(mod, func_name)
            vec_func = np.vectorize(func)    
        else:
            func = getattr(FunctionApprox2d, func_str)
            vec_func = np.vectorize(func)    
        func_values = vec_func(X)
        return func_values

    @staticmethod
    def myfunc(x):
        return x
    @staticmethod
    def x2(x):
        return x**2
    @staticmethod
    def x5(x):
        return x**5
    @staticmethod
    def func1(x):
        return 2/3.14*math.acos(x)
    @staticmethod
    def func3(x):
        return np.abs(math.sin(2*3.14*x))
    @staticmethod
    def func2(x):
        return 1/2*(1+math.sin(2*3.14*x))


    def plot_curve(self):
        plt.plot(self.func_points[:,0], self.func_points[:,1], "g--")
        plt.title(self.func_str)
        for i, edge in enumerate(self.func_path):
            points = self.points[edge]
            plt.plot(points[:,0], points[:,1], "r")
        plt.scatter(self.points[self.nearest_points][:,0], self.points[self.nearest_points][:,1], s=100)
        plt.scatter(self.func_points[:,0], self.func_points[:,1], c="r")
