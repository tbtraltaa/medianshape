import numpy as np 

def interpolate_segment(StartPoint, EndPoint, NumberofInterpolation):
    if len(StartPoint)== 2:
        x = np.linspace(StartPoint[0], EndPoint[0], NumberofInterpolation)
        y = np.linspace(StartPoint[1], EndPoint[1], NumberofInterpolation)
        return np.array([x,y]).T
    elif len(StartPoint)== 3:
        x = np.linspace(StartPoint[0], EndPoint[0], NumberofInterpolation)
        y = np.linspace(StartPoint[1], EndPoint[1], NumberofInterpolation)
        z = np.linspace(StartPoint[2], EndPoint[2], NumberofInterpolation)
        return np.array([x,y,z]).T
