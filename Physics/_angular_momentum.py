import numpy as np

def angular_momenta(radi, velocities, masses):
    result = np.cross(radi, velocities - velocities.mean(axis = 0)) * masses[:, None]
    return result

def angular_momentum(radi, velocities, masses):
    result = np.sum(angular_momenta(radi, velocities, masses), axis = 0)
    return result

def specific_angular_momentum(radi, velocities):
    return angular_momentum(radi, velocities, np.ones(radi.shape[0]))
