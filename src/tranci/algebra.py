import numpy as np

def normalize(zaxis):
    zaxis = np.array(zaxis) # turn to vector
    zaxis = zaxis/np.sqrt(zaxis.dot(zaxis)) # normalize
    return zaxis
