import numpy as np
from matplotlib import pyplot as plt

from DataAccess import ParticleReadConversion_EagleSnapshot, load_catalouge_field, Simulations, SimulationModels, ParticleType, UnitSystem, constants
from Physics import angular_momentum

from assembily_history import physical_centeral_mass_positions, assembily_history
from coordinate_transform import produce_disk_x_basis, transform_coordinates

model = SimulationModels.RECAL
tag = "028_z000p000"

relitive_data_root = ".\\gm_for_mphys"

galaxy_centres = [physical_centeral_mass_positions[sim][constants.tags[-1]] for sim in Simulations]#TODO: check whether cgs or physical code is uncommented!!!

simulations = (Simulations.Early, Simulations.Organic, Simulations.Late)
for i in range(len(simulations)):
    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulations[i], model, relitive_data_root)

    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulations[i], SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulations[i]][tag]["halo"] - 1]
    #selection_radius = r_200
    selection_radius = 9.257*(10**22)# 30 KPc

    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
    net_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    for j in range(len(particle_types)):
        _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_types[j], "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_velocities = snapshot.particle_read_sphere(particle_types[j], "Velocity", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        #particle_velocities -= particle_velocities.mean(axis = 0)
        particle_masses = snapshot.particle_read_sphere(particle_types[j], "Mass", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)

        #net_angular_momentum += np.sum(np.cross(particle_locations_box_adjusted, particle_velocities) * particle_masses[:, None], axis = 0)
        net_angular_momentum += angular_momentum(particle_locations_box_adjusted, particle_velocities, particle_masses)

    rotation_axis_vector = net_angular_momentum / np.linalg.norm(net_angular_momentum, 2)
    rotation_plane_x_vector = produce_disk_x_basis(rotation_axis_vector)

    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.dark_matter)
    particle_type_graph_names = ("Gas", "Star" ,"Dark Matter")
    for j in range(len(particle_types)):
        particle_locations, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_types[j], "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_locations_box_adjusted /= constants.SimulationConstants.get_constants()["CM_PER_MPC"]

        particle_locations_transformed, _ = transform_coordinates(particle_locations_box_adjusted, rotation_plane_x_vector, rotation_axis_vector)

        fig = plt.figure(figsize = (12, 12))

        ax = plt.subplot(2, 2, 1, aspect = "equal")
        ax.set_xlabel("X (MPc)")
        ax.set_ylabel("Y (MPc)")
        ax.scatter(particle_locations_transformed[:, 0], particle_locations_transformed[:, 1], s = 0.01)
        x_lims = ax.get_xlim()
        y_lims = ax.get_ylim()
        axis_lims = (min(x_lims[0], y_lims[0]), max(x_lims[1], y_lims[1]))
        ax.set_xlim(*axis_lims)
        ax.set_ylim(*axis_lims)

        ax = plt.subplot(2, 2, 2, aspect = "equal")
        ax.set_xlabel("Z (MPc)")
        ax.set_ylabel("Y (MPc)")
        ax.scatter(particle_locations_transformed[:, 2], particle_locations_transformed[:, 1], s = 0.01)
        x_lims = ax.get_xlim()
        y_lims = ax.get_ylim()
        axis_lims = (min(x_lims[0], y_lims[0]), max(x_lims[1], y_lims[1]))
        ax.set_xlim(*axis_lims)
        ax.set_ylim(*axis_lims)

        ax = plt.subplot(2, 2, 3, aspect = "equal")
        ax.set_xlabel("X (MPc)")
        ax.set_ylabel("Z (MPc)")
        ax.scatter(particle_locations_transformed[:, 0], particle_locations_transformed[:, 2], s = 0.01)
        x_lims = ax.get_xlim()
        y_lims = ax.get_ylim()
        axis_lims = (min(x_lims[0], y_lims[0]), max(x_lims[1], y_lims[1]))
        ax.set_xlim(*axis_lims)
        ax.set_ylim(*axis_lims)
    
        fig.suptitle("{} Particle Projection of the Assembled Galaxy in the {} Regime".format(particle_type_graph_names[j], Simulations.to_string(simulations[i]).title()))
        plt.show()
