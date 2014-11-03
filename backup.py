from __future__ import absolute_import

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

import distmesh as dm

class Mesh():
    def __init__(self, boundary_box, fixed_points=[], interval_size=None):
        self.fixed_points = fixed_points
        self.edges = []
        self.boundary_box = boundary_box
        if not interval_size:
            self.interval_size = 100

    def set_points(self, fixed_points, boundary_box, interval_size=None):
        self.boundary_box = boundary_box
        self.fixed_points = fixed_points
        if not interval_size:
            self.interval_size = 100
            
    def generate_mesh(self):
        print globals()
        self.points, self.simplices = self.square()
        print self.points, self.simplices
        print "MKT2t", dm.mkt2t(self.simplices)
        self.edges = []
        self.set_edges() 
        #self.points_ordered = self.points[self.points[:,1].argsort()]
        self.sample_size = int(math.sqrt(len(self.points)))
        self.ordered_X = np.sort(np.unique(self.points[:,0]))
        self.sample_step = self.ordered_X.size//self.sample_size
        if self.interval_size > self.sample_step:
            self.interval_size = self.sample_step

    def square(self):
        """Square, with size function point and line sources"""
        print globals()
        dist_function = lambda p: dm.drectangle(p,0,1,0,1)
        return dm.distmesh2d(dist_function, dm.huniform, 0.05, self.boundary_box, self.fixed_points)
        #fd = lambda p: np.sqrt((p**2).sum(1))-1.0
        #return dm.distmesh2d(fd, dm.huniform, 0.06, (-1,-1,1,1))

    def set_edges(self):
        if self.points.shape[1] == 2:
            edges = set()
            for simplex in self.simplices:
                for i in range(len(simplex)):
                    edges.add(tuple(sorted([simplex[i], simplex[(i+1)%len(simplex)]])))
            self.edges = np.array(list(edges))

    def generate_curve(self, func_str=None, interval_size=None):
        if interval_size and interval_size < self.sample_step -1:
            self.interval_size = interval_size
        nearest_points, func_points = self.find_nearest_points(func_str)
        print "Nearest_points:\n", nearest_points
        print "Function points:\n", func_points
        func_path = self.find_path(func_str, nearest_points)
        Mesh.disp_edges(func_path, "Function edges")
        self.plot()
        self.plot_curve(func_str, func_points, nearest_points, func_path)
        path_vector = self.get_path_vector(func_path)
        csr_path = csr_matrix(path_vector)
        print "Path vector:\n", csr_path

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

    def find_path(self, func_str, nearest_points):
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

    def find_neighbors(self, point_idx):
        return self.mesh.vertex_neighbor_vertices[1][self.mesh.vertex_neighbor_vertices[0]\
        [point_idx]:self.mesh.vertex_neighbor_vertices[0][point_idx+1]]

    def orient_simplices(self):
        self.initial_simplicies = self.simplices
        direction = self.right_hand_rule(self.simplices[0])
        simplices = np.array(range(0, len(self.simplices)))
        if direction < 0:
            self.simplices[0] = self.simplices[0,::-1]
        for i, simplex in enumerate(self.simplices):
            simplices = np.delete(simplices, np.where(simplices==i))
            neighbors = self.mesh.neighbors[i]
            for opposit_point in np.where(neighbors >= 0)[0]:
                n_simplex_idx = neighbors[opposit_point]
                if any(simplices==n_simplex_idx):
                    n_simplex = self.simplices[n_simplex_idx]
                    n_boundary = Mesh.boundary(n_simplex)
                    subsimplex= Mesh.boundary(simplex, opposit_point)
                    for n_face in n_boundary:
                        if all((np.array(subsimplex) - np.array(n_face)) == 0):
                            self.simplices[n_simplex_idx] = n_simplex[::-1]
                    simplices = np.delete(simplices, np.where(simplices==n_simplex_idx))
        print self.right_hand_rule(self.simplices[0])

    def orient_simplices_2D(self):
        self.initial_simplicies = self.simplices
        for i, simplex in enumerate(self.simplices):
            if self.right_hand_rule(simplex) < 0: 
                self.simplices[i] = self.simplices[i,::-1]
        print self.right_hand_rule(self.simplices[0])

    def right_hand_rule(self, simplex):
        edges = list()
        edges.append(Mesh.boundary(simplex, 0))
        edges.append(Mesh.boundary(simplex, 1))
        edge_points = self.points[np.array(edges)]
        v1 = edge_points[0][1] - edge_points[0][0]
        v2 = edge_points[1][1] - edge_points[1][0]
        direction = np.cross(v1,v2)
        return direction

    @staticmethod
    def simplex_volume():
        pass
        

    @staticmethod
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
    
    @staticmethod
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
        
    def get_path_vector(self, func_path):
        path_vector = np.zeros(shape=(self.edges.shape[0], 1))
        for i, edge in enumerate(self.edges):
            for path_edge in func_path:
                if all((path_edge - edge) == 0):
                    path_vector[i] = 1
                    break
                elif all((path_edge - np.array(list(reversed(edge)))) == 0):
                    path_vector[i] = -1
                    break
        return path_vector
                
    def to_string(self):
        print "Mesh Points:"
        for i, p in enumerate(self.points):
                print i, p
        print "Point numbers in triangle:"
        for i, t in enumerate(self.simplices):
                print i, t
        Mesh.disp_edges(self.edges, "Mesh edges")

    def plot(self):
        plt.triplot(self.points[:,0], self.points[:,1], self.simplices.copy())
        #plt.scatter(self.points[:,0], self.points[:,1])
        
    def plot_curve(self, func_str, func_points, nearest_points, func_path):
        func_points = np.asarray(func_points)
        plt.plot(func_points[:,0], func_points[:,1], "g--")
        plt.title(func_str)
        for i, edge in enumerate(func_path):
            points = self.points[edge]
            plt.plot(points[:,0], points[:,1], "r")
        plt.scatter(self.points[nearest_points][:,0], self.points[nearest_points][:,1], s=100)
        plt.scatter(func_points[:,0], func_points[:,1], c="r")
        plt.show()

    @staticmethod
    def myfunc(x):
        return x
    @staticmethod
    def x2(x):
        return x**2
    @staticmethod
    def x3(x):
        return x**2



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
    mesh = Mesh((0, 0, 1, 1), [(0,0),(0,1),(1,0),(1,1)]);
    mesh.generate_mesh()
    mesh.to_string()
    functions = ['myfunc', 'x2', 'x3', 'math.atan', 'math.acos']
    functions = ['myfunc']
    for f in functions:
        mesh.generate_curve(f)
    mesh.orient_simplices_2D()
    #mesh.orient_simplices()
