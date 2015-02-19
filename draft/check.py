import numpy as np
from scipy import sparse

if __name__ == "__main__":
    m_edges = 3851
    x1 = sparse.dok_matrix((3851, 1), dtype=np.int8) 
    x2 = np.loadtxt('/home/altaa/dumps_large_mesh_full/x-default-lambda-0.0001.txt')
    x2 = x2[0:m_edges] - x2[m_edges: 2*m_edges]
    x2 = sparse.dok_matrix(x2.reshape(x2.size, 1), dtype=np.int8)
    print x2.shape
    with open('/home/altaa/dumps_large_mesh_full/t_BK.txt', 'r') as f:
        for line in f.readlines():
            data = line.split()
            x1[int(data[0]), 0] = np.int8(data[1])
    print x1.shape
    print x2.toarray()-x1.toarray()



