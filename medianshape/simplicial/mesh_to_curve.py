import matplotlib.pyplot as plt 
import numpy as np
from medianshape.simplicial.gaussian_kernel_convolution import *
from medianshape.simplicial.interpolate_segment import *

def smoothen(sampleCurve):
    interpolateSize = 50
    sigma = 150*0.6
    weight = np.array([0.3,0.4,0.3])
    
    ## intepolate points ## 
    smoothCurve = np.empty((0,sampleCurve.shape[1]))
    for i in range(len(sampleCurve)-1):
        addPoints = interpolate_segment(sampleCurve[i], sampleCurve[i+1], interpolateSize)
        smoothCurve = np.append(smoothCurve, addPoints, axis =0)

    smoothCurve[:,0] = gaussian_kernel_convolution(smoothCurve[:,0], sigma)
    smoothCurve[:,1] = gaussian_kernel_convolution(smoothCurve[:,1], sigma)
    smoothCurve[:,2] = gaussian_kernel_convolution(smoothCurve[:,2], sigma)


    ## Convolute points ## 
    # smoothCurve[:,0] = convolution(smoothCurve[:,0], weight,1)
    # smoothCurve[:,1] = convolution(smoothCurve[:,1], weight,1)

    # ## Chaikins Corner Cutting ##
    # smoothCurve = chaikins_corner_cutting(sampleCurve,3)

    return smoothCurve



