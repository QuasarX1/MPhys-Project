import numpy as np

from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem

from assembily_history import physical_centeral_mass_positions, assembily_history
from data_over_time import produce_simulations_graph

relitive_data_root = ".\\gm_for_mphys"

def func(halo, subhalo, tag, simulation):
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulation][tag]["halo"] - 1]
    #selection_radius = r_200
    selection_radius = 9.257*(10**22)# 30 KPc

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    # Calculate the net angular momentum (axis vector)
    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
    net_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    for j in range(len(particle_types)):
        _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_types[j], "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_velocities = snapshot.particle_read_sphere(particle_types[j], "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        particle_velocities -= particle_velocities.mean(axis = 0)
        particle_masses = snapshot.particle_read_sphere(particle_types[j], "Mass", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)

        net_angular_momentum += np.sum(np.cross(particle_locations_box_adjusted, particle_velocities) * particle_masses[:, None], axis = 0)

    # Make a unit vector
    rotation_axis_vector = net_angular_momentum / np.linalg.norm(net_angular_momentum, 2)

    # Calculate angular momenta and kinetic energy for all stars
    _, star_particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.star, "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    star_particle_velocities = snapshot.particle_read_sphere(ParticleType.star, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    star_particle_velocities[:, 0] -= star_particle_velocities[:, 0].mean()
    star_particle_velocities[:, 1] -= star_particle_velocities[:, 1].mean()
    star_particle_velocities[:, 2] -= star_particle_velocities[:, 2].mean()
    star_particle_masses = snapshot.particle_read_sphere(ParticleType.star, "Mass", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)

    angular_momenta = np.cross(star_particle_locations_box_adjusted, star_particle_velocities) * star_particle_masses[:, None]
    kinetic_energies = 0.5 * star_particle_masses * (star_particle_velocities**2).sum(axis = 1)
    
    kappa_CO = kinetic_energies[angular_momenta.dot(rotation_axis_vector) > 0].sum() / kinetic_energies.sum()

    return kappa_CO

produce_simulations_graph(func, "$\kappa_{CO}$", "Co-Rotating Stellar Kinetic Energy Fraction", log_x = True)

