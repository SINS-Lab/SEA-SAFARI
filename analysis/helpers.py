import numpy as np 
import numpy.matlib as matlib

def rotate(a, b):
    # Ensure are normalized first.
    a = a/np.linalg.norm(a)
    b = b/np.linalg.norm(b)

    v = np.cross(a, b)
    c = a.dot(b)
    c = 1/(1+c)
    I = matlib.eye(3)
    V = matlib.matrix([[0    ,-v[2], v[1]],\
                       [v[2] ,0    ,-v[0]],\
                       [-v[1],v[0] ,0    ]])
    return I + V + V.dot(V) * c

def rotate_all(axis, dir, set):
    # Ensure are normalized first.
    axis = axis/np.linalg.norm(axis)
    dir = dir/np.linalg.norm(dir)
    
