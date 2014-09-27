
word_list = open('/usr/share/dict/words').readlines()
word_list = map(str.strip, word_list)
word_list = [word for word in word_list if len(word) == 3]
word_list = [word for word in word_list if word[0].islower()]
word_list = [word for word in word_list if word.isalpha()]
word_list = map(str.lower, word_list)
print len(word_list)
import numpy as np
word_list = np.asarray(word_list)
print word_list.dtype
word_list.sort()
word_bytes = np.ndarray((word_list.size, word_list.itemsize), dtype='int8', buffer=word_list.data)
print word_bytes
print word_bytes.shape

from scipy.spatial.distance import pdist, squareform
from scipy.sparse import csr_matrix
print word_bytes
hamming_dist = pdist(word_bytes, metric='hamming')
print hamming_dist
print len(hamming_dist)
print 1.5 / word_list.itemsize
print hamming_dist < 1.5 / word_list.itemsize
print squareform(hamming_dist < 1.5 / word_list.itemsize)
graph = csr_matrix(squareform(hamming_dist < 1.5 / word_list.itemsize))
i1 = word_list.searchsorted('ape')
i2 = word_list.searchsorted('man')
print word_list[i1]
print word_list[i2]
from scipy.sparse.csgraph import dijkstra
distances, predecessors = dijkstra(graph, indices=i1, return_predecessors=True)
print distances[i2]
path = []
i = i2
while i != i1:
    path.append(word_list[i])
    i = predecessors[i]
path.append(word_list[i1])
print path[::-1]

