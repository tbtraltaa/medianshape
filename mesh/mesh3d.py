from meshpy.tet import MeshInfo, build

import numpy as np

class Mesh():
    shapes = ['rectangle', 'triangle', 'circle']
    def __init__(self, points=[], facets=[]):
        self.points = points
        self.facets = facets
        self.generate_mesh()

    def set_points(self, shape, x, y, x_range, y_range=None):
        #if shape in shapes:
        pass
            
    def generate_mesh(self):
        self.mesh_info = MeshInfo()
        self.mesh_info.set_points(self.points)
        self.mesh_info.set_facets(self.facets)
        self.mesh = build(self.mesh_info)

    def to_string(self):
        print "Mesh Points:"
        for i, p in enumerate(self.mesh.points):
                print i, p
        print "Point numbers in tetrahedra:"
        for i, t in enumerate(self.mesh.elements):
                print i, t
        self.mesh.write_vtk("test.vtk")
        #self.mesh.dump()

    def list_edges(self):
        print "Point numbers in edges:"
        for edge in self.mesh.edges:
            print edge
        self.mesh.save_edges('edges')

    def generate_curve(self, func=None):
        
        pass
    def plot(self):
        pass

if __name__ == "__main__":
    points = [(0,0,0), (2,0,0), (2,2,0), (0,2,0),
        (0,0,12), (2,0,12), (2,2,12), (0,2,12)] 
    facets = [
        [0,1,2,3],
            [4,5,6,7],
                [0,4,5,1],
                    [1,5,6,2],
                        [2,6,7,3],
                            [3,7,4,0]] 
    mesh = Mesh(points, facets);
    mesh.generate_mesh()
    mesh.to_string()
    mesh.list_edges()
