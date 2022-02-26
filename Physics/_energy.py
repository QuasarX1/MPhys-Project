import numpy as np

from Physics._angular_momentum import angular_momenta

def kinetic_energy(mass, velocity):
    return 0.5 * mass * (velocity**2).sum(axis = 1)

def kappa_CO(particle_radi, particle_velocities, particle_masses):
    particle_angular_momenta = angular_momenta(particle_radi, particle_velocities, particle_masses)
    kinetic_energies = kinetic_energy(particle_masses, particle_velocities)
    corotating_kinetic_energies = 0.5 * particle_angular_momenta[:, 1]**2 / particle_masses / (particle_radi[:, 0]**2 + particle_radi[:, 2]**2)
    return corotating_kinetic_energies[particle_angular_momenta.dot(np.array([0, 1, 0])) > 0].sum() / kinetic_energies.sum()
