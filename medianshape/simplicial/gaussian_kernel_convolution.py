from scipy import ndimage
import numpy as np

def gaussian_kernel_convolution(Data, Sigma):
	return ndimage.gaussian_filter1d(Data, Sigma, mode = 'nearest')

def convolution(Data, Weights, Boundary = 0):
	convoluted = np.convolve(Data, Weights, 'full')
	if Boundary == 0:
		return convoluted
	elif Boundary == 1:
		convoluted[0] = Data[0]
		convoluted[-1] = Data[-1]
		return convoluted


def chaikins_corner_cutting(coords, refinements=2):
    coords = np.array(coords)
    for _ in range(refinements):
        L = coords.repeat(2, axis=0)
        R = np.empty_like(L)
        R[0] = L[0]
        R[2::2] = L[1:-1:2]
        R[1:-1:2] = L[2::2]
        R[-1] = L[-1]
        coords = L * 0.75 + R * 0.25
    return coords