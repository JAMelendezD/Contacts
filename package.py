import numpy as np

def distance(positions,l):
    '''
    Calculates a matrix of distances between vectors with periodic boundaries
    '''
    if type(positions) != type(np.ndarray):
        raise TypeError("Input must be a numpy array")
    if len(positions) == 0:
        raise ValueError("Numpy array is empty")
    
    dists = np.zeros(np.shape(positions))
    for i,pos, in enumerate(positions):
        closest = np.remainder(pos - positions + l/2.0, l) - l/2.0
        distances = np.sqrt(np.einsum("ij,ij->i", closest, closest))
    return distances
