import numpy as np
from numba import jit

def distance_mat(positions,l):
    '''
    Calculates a matrix of distances between vectors with periodic boundaries
    '''
    if type(positions) != type(np.zeros(1)):
        raise TypeError("Input must be a numpy array")
    if type(l) not in [int,float] or l <= 0:
        raise TypeError("Lattice lenght must be a positive real number")
    
    rows, _ = np.shape(positions)
    
    if rows <= 1:
        raise ValueError("Numpy array is empty or contain only one distance")

    distances = np.zeros((rows,rows))
    for i,pos in enumerate(positions):
        closest = np.remainder(pos - positions + l/2.0, l) - l/2.0
        distances[i] = np.sqrt(np.einsum("ij,ij->i", closest, closest))
    return distances


def contacts(positionsA,positionsB,cutoff):
    '''
    Creates a contact map based on a cutoff of the distances
    '''
    subts = positionsA[:,None,:] - positionsB
    distances = np.sqrt(np.einsum('ijk,ijk->ij',subts,subts))
    counter = (distances <= cutoff).astype(int)
    return counter