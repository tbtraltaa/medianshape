import importlib
import random
import math

import numpy as np
import matplotlib.pyplot as plt

from scipy.spatial.distance import pdist

from scipy.sparse.csgraph import shortest_path
from scipy.sparse.csgraph import dijkstra
from scipy.sparse import csr_matrix

from meshpy.triangle import MeshInfo, build, write_gnuplot_mesh


class Mesh():
    shapes = ['rectangle', 'triangle', 'circle']
    def __init__(self, points=[], facets=[], max_volume=100):
        self.points = points
        self.facets = facets
        self.max_volume = max_volume
        self.x = 0
        self.y = 0
        self.x_range = 1
        self.y_range = 1
        self.interval_size = 3

    def set_points(self, x, y, x_range, y_range=None, shape=None):
        self.x = x
        self.y = y
        self.x_range = x_range
        if not y_range:
            self.y_range = x_range
        self.points = [(self.x,self.y), (self.x, self.y_range), (self.x_range, self.y ), (self.x_range,self.y_range)]
        self.facets = [(0,1), (0,2), (1,3), (2,3)]
            
    def generate_mesh(self):
        self.mesh_info = MeshInfo()
        self.mesh_info.set_points(self.points)
        self.mesh_info.set_facets(self.facets)
        self.mesh = build(self.mesh_info, generate_faces=True, min_angle = 30, max_volume=self.max_volume)

        self.mesh_points = np.array(self.mesh.points).reshape(len(self.mesh.points), 2)
        mesh_points_idx = np.array(range(0, len(self.mesh.points))).reshape(len(self.mesh.points), 1)
        self.mesh_points = np.concatenate((mesh_points_idx, self.mesh_points), axis=1)

        self.mesh_faces = np.array(self.mesh.faces).reshape(len(self.mesh.faces), 2)
        mesh_faces_idx = np.array(range(0, len(self.mesh.faces))).reshape(len(self.mesh.faces), 1)
        self.mesh_faces = np.hstack((mesh_faces_idx, self.mesh_faces))

        self.mesh_elements = np.array(self.mesh.elements).reshape(len(self.mesh.elements), 3)
        mesh_elements_idx = np.array(range(0, len(self.mesh.elements))).reshape(len(self.mesh.elements), 1)
        self.mesh_elements = np.hstack((mesh_elements_idx, self.mesh_elements))
        self.point_X= np.sort(np.unique(self.mesh_points[:,1]))

        #self.mesh_points_ordered = self.mesh_points[self.mesh_points[:,1].argsort()]

    def to_string(self):
        print "Mesh Points:"
        for i, p in enumerate(self.mesh.points):
                print i, p
        print "Point numbers in triangle:"
        for i, t in enumerate(self.mesh.elements):
                print i, t

    def list_edges(self):
        print "Point numbers in edges:"
        for i, face in enumerate(self.mesh.faces):
            print i, face
        write_gnuplot_mesh("faces", self.mesh)

    def generate_curve(self, func_str=None, interval_size=None):
        if interval_size:
            self.interval_size = interval_size
        sample_points = random.sample(self.mesh_points, int(math.sqrt(self.mesh_points.size)))
        sample_points = np.array(sample_points)
        sample_X = np.unique(sample_points[:, 1])
        sample_X = sample_X.reshape(sample_X.size, 1)
        print "Sample points"
        print sample_points
        optimum_points= self.find_closest_points(func_str, sample_X)
        print "Optimum_points\n", optimum_points
        func_values = Mesh.vectorize_func(func_str, sample_X) 
        func_values = func_values.reshape(func_values.size, 1)
        func_points = np.concatenate((sample_X, func_values), axis=1)
        func_edges = []
        print "Function points:"
        print func_points
        func_edges = self.find_func_edges(func_str, optimum_points)
        Mesh.disp_edges(func_edges, "Function edges")
        self.plot()
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

    @staticmethod
    def find_diff(func_str, points):
        metric = abs(sum(points[:,2] - Mesh.vectorize_func(func_str, points[:,1])))
        #metric += pdist(points[:,1:], 'euclidean')
        return metric

    def find_closest_points(self, func_str, sample_X):
        sample_X = np.sort(sample_X)
        optimum_points = list()
        for i, x in enumerate(sample_X):
            optimum_points.append(self.find_optimum_point(func_str, x))
        return  np.array(optimum_points)

    def find_optimum_point(self, func_str, x):
        candidate_points = self.mesh_points[np.where(self.mesh_points[:,1]==x)]        
        interval_X = self.find_interval_X(x)
        interval_func_values = Mesh.vectorize_func(func_str, interval_X)
        min_diff = 1000
        optimum_point = list()
        for point in candidate_points:
            d = interval_func_values - point[2]
            diff = abs(sum(interval_func_values - point[2])) 
            if diff < min_diff:
                min_diff = diff
                optimum_point = point.tolist()
        return optimum_point

    def find_interval_X(self, x):
        x_idx = np.where(self.point_X==x)[0]
        interval_X = []      
        i = 0
        start_idx = 0
        while True:
            interval_idx = x_idx - self.interval_size//2 + i
            if interval_idx >= 0:
                start_idx = interval_idx
                break
            i += 1
        interval_X = self.point_X[start_idx:start_idx + self.interval_size]
        return interval_X

    def find_closest_edge(self, prev_edge_idx, curr_point, next_point=[]):
        edges = {}
        closest_edge_idx = -1
        min_dist = -1
        first_edge = True
        if prev_edge_idx !=None:
            edges = self.find_connected_edges(prev_edge_idx)
        else:
            for i, edge in enumerate(self.mesh.faces):
                edges[i] = edge
        
        for i, edge in edges.iteritems():
            dist = Mesh.dist(self.mesh.points[edge[0]], self.mesh.points[edge[1]], curr_point) 
            if len(next_point) != 0 :
                dist += Mesh.dist(self.mesh.points[edge[0]], self.mesh.points[edge[1]], next_point) 
            if first_edge:
                min_dist = dist
                closest_edge_idx = i
                first_edge = False
            if dist < min_dist:
                min_dist = dist
                closest_edge_idx = i
        return  closest_edge_idx

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
        print type(self.mesh_points[:,0])
        print type(self.mesh_points[:,1])
        print type(self.mesh_elements.copy())
        for i, edge in enumerate(self.mesh.faces):
            X = []
            Y = []
            for point in edge:
                X.append(self.mesh.points[point][0])
                Y.append(self.mesh.points[point][1])
        plt.plot(X, Y)
        plt.plot(self.mesh_points[:,1], self.mesh_points[:,2], 'ro')
        #plt.triplot(self.mesh_points[:,1], self.mesh_points[:,2], self.mesh_elements)


    def plot_curve(self, func_points, func_edges, optimum_points):
        func_points = np.asarray(func_points)
        plt.plot(func_points[:,0], func_points[:,1], "g--")
        X = []
        Y = []
        for i, edge in enumerate(func_edges):
            X = []
            Y = []
            for point in edge:
                X.append(self.mesh.points[point][0])
                Y.append(self.mesh.points[point][1])
            plt.plot(X, Y, "r")
        plt.scatter(optimum_points[:,1], optimum_points[:,2], s=100)
        plt.show()

if __name__ == "__main__":
    #points = [(0,0), (10, 0), (10, 10), (0, 10)]
    #facets = [(0,1), (1,2), (2,3), (3,0)]
    mesh = Mesh();
    mesh.set_points(0, 0, 50)
    mesh.generate_mesh()
    mesh.to_string()
    mesh.list_edges()
    mesh.generate_curve("myfunc")
