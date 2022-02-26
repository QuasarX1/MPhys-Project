import numpy as np

def angular_momenta(radi, velocities, masses):
    return np.cross(radi, velocities - velocities.mean(axis = 0)) * masses[:, None]

def angular_momentum(radi, velocities, masses):
    return np.sum(angular_momenta(radi, velocities, masses), axis = 0)

def specific_angular_momentum(radi, velocities):
    return angular_momentum(radi, velocities, np.ones(radi.shape[0]))
