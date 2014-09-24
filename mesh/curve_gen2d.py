import importlib
import random
import math

from sets import Set

import numpy as np
import numpy.ma as ma
from meshpy.triangle import MeshInfo, build, write_gnuplot_mesh

import matplotlib.pyplot as plt

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
        self.mesh = build(self.mesh_info, generate_faces=True, max_volume=self.max_volume)
        print type(self.mesh.points)
        self.mesh_points = np.array(np.array(self.mesh.points), dtype=[("x", float), ("y", float)]).reshape(len(self.mesh.points), 2)
        self.mesh_point_indices = np.array(np.array(range(0, len(self.mesh.points))), dtype=[("idx", float)]).reshape(len(self.mesh.points), 1)
        print type(self.mesh_point_indices), self.mesh_point_indices.shape
        print type(self.mesh_points), self.mesh_points.shape

        print self.mesh_points.x
        self.mesh_points = np.concatenate((self.mesh_point_indices, self.mesh_points), axis=1)

        self.mesh_faces = np.array(self.mesh.faces).reshape(len(self.mesh.faces), 2)
        self.mesh_face_indices = np.array(range(0, len(self.mesh.faces))).reshape(len(self.mesh.faces), 1)
        self.mesh_faces = np.hstack((self.mesh_face_indices, self.mesh_faces))

        self.mesh_elements = np.array(self.mesh.elements).reshape(len(self.mesh.elements), 3)
        self.mesh_element_indices = np.array(range(0, len(self.mesh.elements))).reshape(len(self.mesh.elements), 1)
        self.mesh_elements = np.hstack((self.mesh_element_indices, self.mesh_elements))
        self.mesh_points_ordered = np.sort(self.mesh_points, order="x")
        print self.mesh_points_ordered


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

        for points in self.mesh.elements:
            for pt in points:
                pass
                #print "%f %f" % tuple(self.mesh.points[pt])
            #print "\n"

    @staticmethod
    def disp_edges(edges):
        print "Point numbers in edges:"
        for i, edge in enumerate(edges):
            print i, edge
    
    def generate_curve(self, func_str=None):
        sample_points = random.sample(self.mesh_points, int(math.sqrt(self.mesh_points.size)))
        sample_points = np.array(sample_points)
        print "Sample points"
        print sample_points
        #sample_points = np.concatenate((sample_Xs, self.find_closest_point(func_str, sample_Xs)), axis=1)
        self.find_closest_point(func_str, sample_points)

        func_points = []
        func_edges = []
        print "Function points:"
        print func_points
        Mesh.disp_edges(func_edges)
        self.plot()
        self.plot_curve(func_points, func_edges)

    def find_closest_point(self, func_str, sample_points):
        sample_Xs = np.sort(np.unique(sample_points[:,1]))
        func_values = Mesh.vectorize_func(func_str, sample_Xs)[:,1] 
        print sample_Xs
        for i, x in enumerate(sample_Xs):
            y = self.mesh_points[self.mesh_points[:,1]==x][:,2]
            #print  abs(func_values[i] - y)
        return []

    def find_interval_points(self, x):
        pass
         
        

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

    def find_connected_edges(self, edge_idx):
        connected_edges = {}
        for i, edge in enumerate(self.mesh.faces):
            if self.mesh.faces[edge_idx][1] == edge[0]:
                connected_edges[i]= edge
        return connected_edges

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
        Y = vec_func(X)
        X = X.reshape(X.size, 1)
        Y = Y.reshape(Y.size, 1)
        func_points = np.concatenate((X, Y), axis=1)
        return func_points

    def plot(self):
        plt.figure(0)
        X = []
        Y = []
        for i, edge in enumerate(self.mesh.faces):
            for point in edge:
                X.append(self.mesh.points[point][0])
                Y.append(self.mesh.points[point][1])
            plt.plot(X, Y)

        X = []
        Y = []
        for i, p in enumerate(self.mesh.points):
            X.append(p[1])
            Y.append(p[0])
        plt.scatter(X,Y)
        plt.show()

    def plot_curve(self, func_points, func_edges):
        func_points = np.asarray(func_points)
        plt.scatter(func_points[:,0], func_points[:,1], c="r")
        for i, edge in enumerate(func_edges):
            X = []
            Y = []
            for point in edge:
                X.append(self.mesh.points[point][0])
                Y.append(self.mesh.points[point][1])
            plt.plot(X, Y, "r--")
        plt.show()

def ufunction(x):
    return x
    
if __name__ == "__main__":
    #points = [(0,0), (10, 0), (10, 10), (0, 10)]
    #facets = [(0,1), (1,2), (2,3), (3,0)]
    mesh = Mesh();
    mesh.set_points(0, 0, 100)
    mesh.generate_mesh()
    mesh.to_string()
    mesh.list_edges()
    mesh.generate_curve("myfunc")
