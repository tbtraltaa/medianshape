import importlib

import numpy as np
from meshpy.triangle import MeshInfo, build, write_gnuplot_mesh

import matplotlib.pyplot as plt

class Mesh():
    shapes = ['rectangle', 'triangle', 'circle']
    def __init__(self, points=[], facets=[], max_volume=50):
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

    def to_string(self):
        print "Mesh Points:"
        for i, p in enumerate(self.mesh.points):
                print i, p
        print "Point numbers in triangle:"
        for i, t in enumerate(self.mesh.elements):
                print i, t
        #self.mesh.write_vtk("test.vtk")
        #self.mesh.dump()

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
        if func_str:
            func_points = Mesh.vectorize_func(func_str, self.x, self.x_range)
            func_edges = []
            edge_idx = None
            for i, point in enumerate(func_points):
                if i+1 < len(func_points):
                    edge_idx = self.find_closest_edge(edge_idx, point, func_points[i+1])
                    func_edges.append(self.mesh.faces[edge_idx])
                else:
                    edge_idx = self.find_closest_edge(edge_idx, point)
                    func_edges.append(self.mesh.faces[edge_idx])
            print "Function points:"
            print func_points
            Mesh.disp_edges(func_edges)
            self.plot()
            self.plot_curve(func_points, func_edges)

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
    def vectorize_func(func_str, x, x_range):
        func_points = []
        if func_str.find(".") != -1:
            mod_name, func_name = func_str.rsplit('.', 1)
            mod = importlib.import_module(mod_name)
            func = getattr(mod, func_name)
            vec_func = np.vectorize(func)    
        else:
            func = getattr(Mesh, func_str)
            vec_func = np.vectorize(func)    
        X = np.array(range(x, x + x_range, 1))
        Y = vec_func(X)
        X = X.reshape(X.size, 1)
        Y = Y.reshape(Y.size, 1)
        func_points = np.concatenate((X, Y), axis=1)
        return func_points

    # Returns distance between p3 and the line containing p1 and p2
    @staticmethod
    def dist(p1, p2, p3):
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
        dividend = np.linalg.norm(np.cross((p2 - p1),(p1 - p3)))
        divisor = np.linalg.norm(p2-p1)
        dist =  dividend/divisor
        return dist 

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