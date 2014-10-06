import importlib
import random
import math

import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from scipy.spatial.distance import pdist, cdist
import matplotlib.pyplot as plt

from scipy.sparse.csgraph import shortest_path
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix

class Mesh():
    def __init__(self, points=[], x_range=1, y_range=None, interval_size=None):
        self.points = points
        self.x_range = x_range
        self.y_range = y_range
        if not self.y_range:
            self.y_range = x_range
        if not interval_size:
            self.interval_size = 3

    def set_points(self, points, x_range=1, y_range=None, interval_size=None):
        self.x_range = x_range
        self.y_range = y_range
        if not self.y_range:
            self.y_range = x_range
        self.points = points
        if not interval_size:
            self.interval_size = 3
            
    def generate_mesh(self):
        self.mesh = Delaunay(self.points)
        self.edges = []
        self.set_edges() 
        #self.mesh.points_ordered = self.mesh.points[self.mesh.points[:,1].argsort()]
        self.sample_size = int(math.sqrt(len(self.mesh.points)))
        self.ordered_X = np.sort(np.unique(self.mesh.points[:,0]))
        self.sample_step = self.ordered_X.size//self.sample_size
        if self.interval_size > self.sample_step:
            self.interval_size = self.sample_step - 1

    def set_edges(self):
        if self.mesh.points.shape[1] == 2:
            for i, simplex in enumerate(self.mesh.simplices):
                edge0 = [simplex[0], simplex[1]]
                edge1 = [simplex[1], simplex[2]]
                edge2 = [simplex[2], simplex[0]]
                if not self.is_in_edges(edge0):
                    self.edges.append(edge0) 
                if not self.is_in_edges(edge1):
                    self.edges.append(edge1) 
                if not self.is_in_edges(edge2):
                    self.edges.append(edge2) 
        self.edges = np.array(self.edges)

    def is_in_edges(self, edge):
        for e in self.edges:
            if (e[0] == edge[0] and e[1] == edge[1]) or (e[0] == edge[1] and e[1] == edge[0]):
                return True
        return False

    def generate_curve(self, func_str=None, interval_size=None):
        if interval_size and interval_size < self.sample_step -1:
            self.interval_size = interval_size
        #print "Sample points:\n", sample_points
        nearest_points, func_points = self.find_nearest_points(func_str)
        print "Nearest_points:\n", nearest_points
        func_edges = []
        print "Function points:\n", func_points
        func_edges = self.find_path(func_str, nearest_points)
        Mesh.disp_edges(func_edges, "Function edges")
        self.plot()
        self.plot_curve(func_points, nearest_points, func_edges)

    def find_path(self, func_str, nearest_points):
        no_of_points = self.mesh.points.shape[0] 
        adjacency_matrix = np.zeros((no_of_points, no_of_points), dtype=np.float64)  
        for i, p1 in enumerate(self.mesh.points):
            neighbor_indices = self.find_neighbors(i) 
            for j  in neighbor_indices:
                p1 = p1.T.reshape(1,2)
                p2 = self.mesh.points[j].T.reshape(1,2)
                adjacency_matrix[i,j] = 1 + cdist(p1, p2, 'euclidean')[0]

        #graph = csr_matrix(adjacency_matrix)  
        graph = adjacency_matrix
        #print "Adjacent matrix \n", adjacency_matrix
        #print graph.data
        path = list()
        for i, point in reversed(list(enumerate(nearest_points))): 
            print i, point
            if i-1 >= 0:
                i1 = nearest_points[i-1]
                i2 = point
                distances, predecessors = dijkstra(graph, directed=False, \
                                            indices=i1, return_predecessors=True)
                j = i2
                while j != i1:
                    if len(path) == 0:
                        path.append(j)
                        graph[i][j] = 100
                        #graph[j][1] = graph[j][0][1] + 100
                        #print graph[j]

                    elif len(path) > 0 and path[-1]!= j:
                        path.append(j)
                        graph[i][j] = 100
                        #print graph[j]
                        #graph[j][1] = graph[j][0][1] + 100
                    j = predecessors[j]
                path.append(i1)
        faces = list()
        path = np.array(path)
        print "Path", path
        for i, p in reversed(list(enumerate(path))):
            if i+1 < len(path):
                tmp = self.edges == (p, path[i+1])
                edge1 = np.any(np.logical_and(tmp[:,0], tmp[:,1]))
                if edge1:
                    faces.append([int(p), int(path[i+1])])
                tmp = self.edges == (path[i+1], p)
                edge2 = np.any(np.logical_and(tmp[:,0], tmp[:,1]))
                if edge2:
                    faces.append([int(path[i+1]), int(p)])
        return  faces

    def find_neighbors(self, point_idx):
        return self.mesh.vertex_neighbor_vertices[1][self.mesh.vertex_neighbor_vertices[0]\
        [point_idx]:self.mesh.vertex_neighbor_vertices[0][point_idx+1]]

    def find_nearest_points(self, func_str):
        sample_X = []
        nearest_points = list()
        for i in range(0, self.ordered_X.size, self.sample_step):
            nearest_points.append(self.find_nearest_point(func_str, i))
            sample_X.append(self.ordered_X[i])
        func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), \
        Mesh.vectorize(func_str, sample_X).reshape(len(sample_X),1)))
        return  np.array(nearest_points), func_points

    def find_nearest_point(self, func_str, x_idx):
        interval_X = self.find_interval_X(x_idx)
        func_values = Mesh.vectorize(func_str, interval_X)
        min_dist = 1000
        nearest_point = -1
        func_points = np.hstack((interval_X, func_values))
        for i, x in enumerate(interval_X):
            #Assigning candidate indices out of a tuple with 2 elements where the second one is empty
            candidate_idx = np.where(self.mesh.points[:,0]==x)[0]     
            candidate_points = self.mesh.points[candidate_idx]
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

    def to_string(self):
        print "Mesh Points:"
        for i, p in enumerate(self.mesh.points):
                print i, p
        print "Point numbers in triangle:"
        for i, t in enumerate(self.mesh.simplices):
                print i, t
        Mesh.disp_edges(self.edges, "Mesh edges")

    def plot(self):
        plt.triplot(self.mesh.points[:,0], self.mesh.points[:,1], self.mesh.simplices.copy())
        plt.plot(self.mesh.points[:,0], self.mesh.points[:,1], 'yo')
        
    def plot_curve(self, func_points, nearest_points, func_edges):
        func_points = np.asarray(func_points)
        plt.plot(func_points[:,0], func_points[:,1], "g--")
        for i, edge in enumerate(func_edges):
            points = self.mesh.points[edge]
            plt.plot(points[:,0], points[:,1], "r")
        plt.scatter(self.mesh.points[nearest_points][:,0], self.mesh.points[nearest_points][:,1], s=100)
        plt.scatter(func_points[:,0], func_points[:,1], c="r")
        plt.show()

    @staticmethod
    def myfunc(x):
        return x*x

    @staticmethod
    def vectorize(func_str, X):
        func_points = []
        if func_str.find(".") != -1:
            mod_name, func_name = func_str.rsplit('.', 1)
            mod = importlib.import_module(mod_name)
            func = getattr(mod, func_name)
            vec_func = np.vectorize(func)    
        else:
            func = getattr(Mesh, func_str)
            vec_func = np.vectorize(func)    
        func_values = vec_func(X)
        return func_values
    
    @staticmethod
    def disp_edges(edges, title):
        print title
        for i, edge in enumerate(edges):
            print i, edge

if __name__ == "__main__":
    mesh = Mesh();
    #mesh.set_points(np.array([[0,0],[0,1], [1,1], [1,0], [0.5, 0.5],[0.7, 0.7],[0.9, 0.9],[0.3, 0.3]]))
    #mesh.set_points(np.append(np.random.uniform(size=(50,2)),[[0,0],[0,1],[1,0],[1,1]], axis=0))
    mesh.set_points(np.append(np.random.rand(25,2),[[0,0],[0,1],[1,0],[1,1]], axis=0))
    mesh.generate_mesh()
    mesh.to_string()
    mesh.generate_curve("myfunc")
