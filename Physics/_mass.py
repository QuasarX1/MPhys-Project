import numpy as np

def centre_of_mass(positions: np.ndarray, masses: np.ndarray):
    single_result = False
    if len(positions.shape) < 2:
        positions = positions.reshape(positions.shape[0], 1)
    com = np.sum(positions * masses[:, None], axis = 0) / sum(masses)
    return com[0] if single_result else com
