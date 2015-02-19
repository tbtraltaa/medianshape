import numpy as np
points = np.array([[0, 0], [0, 1.1], [1, 0], [1, 1]])
#points = np.array([[0,0,0], [0,1,0], [1,0,0], [0,0,1]])
from scipy.spatial import Delaunay
tri = Delaunay(points)

import matplotlib.pyplot as plt
tri = Delaunay(points)

import matplotlib.pyplot as plt
plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
plt.show()

for i, point in enumerate(tri.points):
    print i, point
print "Simplices\n", tri.simplices
print "Simplices points\n", points[tri.simplices]

print "Neighbors of 1\n", tri.neighbors[1]
print "simplices 1,1\n", tri.simplices[1,1]
print points[tri.simplices[1,1]]
#p = np.array([(0.1, 0.2), (1.5, 0.5)])
p = np.array([(0, 0), (1, 0), (0,1)])

print "p", p

print "Simplices containing P\n", tri.find_simplex(p)
