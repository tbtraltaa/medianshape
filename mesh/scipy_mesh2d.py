import importlib
import random
import math

import numpy as np
from scipy.spatial import Delaunay, ConvexHull
from scipy.spatial.distance import pdist
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
            self.interval_size = 4

    def set_points(self, points, x_range=1, y_range=None, interval_size=None):
        self.x_range = x_range
        self.y_range = y_range
        if not self.y_range:
            self.y_range = x_range
        #points = np.append(points,[(0,0), (0, self.y_range), (self.x_range, 0 ), (self.x_range,0)], axis=0)
        self.points = points
        if not interval_size:
            self.interval_size = 4
            
    def generate_mesh(self):
        self.mesh = Delaunay(self.points)
        self.edges = []
        self.mesh_points = np.array(self.mesh.points).reshape(len(self.mesh.points), 2)
        mesh_points_idx = np.array(range(0, len(self.mesh.points))).reshape(len(self.mesh.points), 1)
        self.mesh_points = np.concatenate((mesh_points_idx, self.mesh_points), axis=1)

        self.set_edges() 

        self.mesh_faces = np.array(self.edges).reshape(len(self.edges), 2)
        mesh_faces_idx = np.array(range(0, len(self.edges))).reshape(len(self.edges), 1)
        self.mesh_faces = np.hstack((mesh_faces_idx, self.edges))

        self.mesh_elements = np.array(self.mesh.simplices).reshape(len(self.mesh.simplices), 3)
        mesh_elements_idx = np.array(range(0, len(self.mesh.simplices))).reshape(len(self.mesh.simplices), 1)
        self.mesh_elements = np.hstack((mesh_elements_idx, self.mesh_elements))
        #self.mesh_points_ordered = self.mesh_points[self.mesh_points[:,1].argsort()]
        self.sample_size = int(math.sqrt(len(self.mesh_points)))
        self.ordered_X = np.sort(np.unique(self.mesh_points[:,1]))
        self.sample_step = self.ordered_X.size//self.sample_size
        if self.interval_size > self.sample_step:
            self.interval_size = self.sample_step

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

    def to_string(self):
        print "Mesh Points:"
        for i, p in enumerate(self.mesh.points):
                print i, p
        print "Point numbers in triangle:"
        for i, t in enumerate(self.mesh.simplices):
                print i, t
        Mesh.disp_edges(self.edges, "Mesh edges")

    def generate_curve(self, func_str=None, interval_size=None):
        if interval_size:
            self.interval_size = interval_size
        #print "Sample points:\n", sample_points
        optimum_points, func_points = self.find_closest_points(func_str)
        print "Optimum_points:\n", optimum_points
        func_edges = []
        print "Function points:\n",func_points
        func_edges = self.find_func_edges(func_str, optimum_points)
        Mesh.disp_edges(func_edges, "Function edges")
        self.plot()
        plt.scatter(func_points[:,0], func_points[:,1], c="r")
        plt.scatter(optimum_points[:,1], optimum_points[:,2], c="y")
        self.plot_curve(func_points, func_edges, optimum_points)

    def find_func_edges(self, func_str, optimum_points):
        no_of_points = self.mesh_points.shape[0] 
        adjacency_matrix = np.zeros((no_of_points, no_of_points), dtype=int)  
        for i, p1 in enumerate(self.mesh_points):
            for j, p2 in enumerate(self.mesh_points): 
                tmp = self.mesh_faces[:, 1:] == (p1[0], p2[0])
                edge1 = np.any(np.logical_and(tmp[:,0], tmp[:,1]))
                tmp = self.mesh_faces[:, 1:] == (p2[0], p1[0])
                edge2 = np.any(np.logical_and(tmp[:,0], tmp[:,1]))
                if edge1 or edge2:
                    p1 = p1.reshape(p1.size,1)
                    p2 = p2.reshape(p2.size,1)
                    adjacency_matrix[i,j] = 1 + Mesh.find_diff(func_str, np.append(p1.T, p2.T, axis=0))

        graph = csr_matrix(adjacency_matrix)  
        print "Adjacent matrix \n", adjacency_matrix
        path = list()
        for i, point in enumerate(optimum_points): 
            if i+1 < optimum_points.shape[0]:
                i1 = point[0]
                i2 = optimum_points[i+1][0]
                distances, predecessors = dijkstra(graph, indices=i1, return_predecessors=True)
                j = i2
                while j != i1:
                    path.append(self.mesh_points[j].tolist())
                    j = predecessors[j]
                path.append(self.mesh_points[i1].tolist())
        faces = list()
        for i, p in reversed(list(enumerate(path))):
            if i+1 < len(path):
                tmp = self.mesh_faces[:, 1:] == (p[0], path[i+1][0])
                edge1 = np.any(np.logical_and(tmp[:,0], tmp[:,1]))
                if edge1:
                    faces.append([int(p[0]), int(path[i+1][0])])
                tmp = self.mesh_faces[:, 1:] == (path[i+1][0], p[0])
                edge2 = np.any(np.logical_and(tmp[:,0], tmp[:,1]))
                if edge2:
                    faces.append([int(path[i+1][0]), int(p[0])])
        return  faces

    def find_closest_points(self, func_str):
        sample_X = []
        optimum_points = list()
        for i in range(0, self.ordered_X.size, self.sample_step):
            optimum_points.append(self.find_optimum_point(func_str, i))
            sample_X.append(self.ordered_X[i])
        func_points = np.hstack((np.array(sample_X).reshape(len(sample_X),1), Mesh.vectorize_func(func_str, sample_X).reshape(len(sample_X),1)))
        return  np.array(optimum_points), func_points

    def find_optimum_point(self, func_str, x_idx):
        interval_X = self.find_interval_X(x_idx)
        func_values = Mesh.vectorize_func(func_str, interval_X)
        min_diff = 1000
        optimum_point = list()
        func_points = np.hstack((interval_X, func_values))
        for i, x in enumerate(interval_X):
            candidate_points = self.mesh_points[np.where(self.mesh_points[:,1]==x)]        
            for cpoint in candidate_points:
                points = np.vstack((func_points[i], cpoint[-2:])) 
                diff = pdist(points) 
                if diff < min_diff:
                    min_diff = diff
                    optimum_point = cpoint.tolist()
        return optimum_point

    @staticmethod
    def find_diff(func_str, points):
        metric = pdist(points[:,-2:], 'euclidean')
        return metric

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
        interval_X = self.ordered_X[start_idx:start_idx + self.interval_size]
        return interval_X.reshape(interval_X.size, 1)

    @staticmethod
    def myfunc(x):
        return x

    @staticmethod
    def vectorize_func(func_str, X):
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

    def plot(self):
        plt.triplot(self.mesh.points[:,0], self.mesh.points[:,1], self.mesh.simplices.copy())
        plt.plot(self.mesh.points[:,0], self.mesh.points[:,1], 'o')
        
    def plot_curve(self, func_points, func_edges, optimum_points):
        func_points = np.asarray(func_points)
        plt.plot(func_points[:,0], func_points[:,1], "g--")
        X = []
        Y = []
        for i, edge in enumerate(func_edges):
            for point in edge:
                X.append(self.mesh.points[point][0])
                Y.append(self.mesh.points[point][1])
            plt.plot(X, Y, "r")
        plt.scatter(optimum_points[:,1], optimum_points[:,2], s=100)
        plt.show()

if __name__ == "__main__":
    mesh = Mesh();
    #mesh.set_points(np.array([[0,0],[0,1], [1,1], [1,0], [0.5, 0.5],[0.7, 0.7],[0.9, 0.9],[0.3, 0.3]]))
    mesh.set_points(np.append(np.random.uniform(size=(12,2)),[[0,0],[0,1],[1,0],[1,1]], axis=0))
    mesh.generate_mesh()
    mesh.to_string()
    mesh.generate_curve("myfunc")
