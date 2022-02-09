import numpy as np
from matplotlib import pyplot as plt
import os
import glob
import cv2

from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, Simulations, SimulationModels, ParticleType, UnitSystem

from assembily_history import physical_centeral_mass_positions, assembily_history
from coordinate_transform import produce_disk_x_basis, transform_coordinates

model = SimulationModels.RECAL

relitive_data_root = ".\\gm_for_mphys"

image_folder = ".\\image_dump"
if not os.path.exists(image_folder):
    os.mkdir(image_folder)



particle_types = (ParticleType.gas, ParticleType.star, ParticleType.dark_matter)
particle_type_graph_names = ("Gas", "Star" ,"Dark Matter")

simulations = (Simulations.Early, Simulations.Organic, Simulations.Late)
for i, simulation in enumerate(simulations):
    for tag in constants.tags:
        snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, model, relitive_data_root)

        potential_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!

        # Skip this tag if there isn't the data present
        if potential_centre is None or assembily_history[simulation][tag]["halo"] is None:
            continue

        r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)[assembily_history[simulation][tag]["halo"] - 1]
        #selection_radius = r_200
        #selection_radius = r_200 / 8#TODO: why is the R_200 too large???
        selection_radius = 9.257*(10**22)# 30 KPc

        axis_calculation_particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
        net_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
        axis_calculation_selection_radius = 9.257*(10**22)# 30 KPc
        for particle_type in axis_calculation_particle_types:
            try:
                _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "Coordinates", potential_centre, axis_calculation_selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
                particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", potential_centre, axis_calculation_selection_radius, UnitSystem.cgs, UnitSystem.cgs)
                particle_velocities -= particle_velocities.mean(axis = 0)
                particle_masses = snapshot.particle_read_sphere(particle_type, "Mass", potential_centre, axis_calculation_selection_radius, UnitSystem.cgs, UnitSystem.cgs)
            except KeyError:
                continue# No particles of this type avalible

            net_angular_momentum += np.sum(np.cross(particle_locations_box_adjusted, particle_velocities) * particle_masses[:, None], axis = 0)

        rotation_axis_vector = net_angular_momentum / np.linalg.norm(net_angular_momentum, 2)
        rotation_plane_x_vector = produce_disk_x_basis(rotation_axis_vector)

        for j in range(len(particle_types)):
            try:
                _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_types[j], "Coordinates", potential_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
            except KeyError:
                continue# No particles of this type avalible
            particle_locations_box_adjusted /= constants.SimulationConstants.get_constants()["CM_PER_MPC"]

            particle_locations_transformed, _ = transform_coordinates(particle_locations_box_adjusted, rotation_plane_x_vector, rotation_axis_vector)

            fig = plt.figure(figsize = (12, 12))

            ax = plt.subplot(2, 2, 1, aspect = "equal")
            ax.set_xlabel("X (pMpc)")
            ax.set_ylabel("Y (pMpc)")
            ax.scatter(particle_locations_transformed[:, 0], particle_locations_transformed[:, 1], s = 0.01)
            x_lims = ax.get_xlim()
            y_lims = ax.get_ylim()
            axis_lims = (min(x_lims[0], y_lims[0]), max(x_lims[1], y_lims[1]))
            ax.set_xlim(*axis_lims)
            ax.set_ylim(*axis_lims)

            ax = plt.subplot(2, 2, 2, aspect = "equal")
            ax.set_xlabel("Z (pMpc)")
            ax.set_ylabel("Y (pMpc)")
            ax.scatter(particle_locations_transformed[:, 2], particle_locations_transformed[:, 1], s = 0.01)
            x_lims = ax.get_xlim()
            y_lims = ax.get_ylim()
            axis_lims = (min(x_lims[0], y_lims[0]), max(x_lims[1], y_lims[1]))
            ax.set_xlim(*axis_lims)
            ax.set_ylim(*axis_lims)

            ax = plt.subplot(2, 2, 3, aspect = "equal")
            ax.set_xlabel("X (pMpc)")
            ax.set_ylabel("Z (pMpc)")
            ax.scatter(particle_locations_transformed[:, 0], particle_locations_transformed[:, 2], s = 0.01)
            x_lims = ax.get_xlim()
            y_lims = ax.get_ylim()
            axis_lims = (min(x_lims[0], y_lims[0]), max(x_lims[1], y_lims[1]))
            ax.set_xlim(*axis_lims)
            ax.set_ylim(*axis_lims)
    
            fig.suptitle("{} Particle Projection of the Assembled Galaxy in the {} Regime".format(particle_type_graph_names[j], Simulations.to_string(simulation).title()))
            plt.savefig(os.path.join(image_folder, "{} {} {}.png".format(Simulations.to_string(simulation).title(), particle_type_graph_names[j], tag)))

    # Create Video
    frameSize = (1200, 1200)
    for j in range(len(particle_types)):
        out = cv2.VideoWriter(os.path.join(image_folder, "{} {}.avi".format(Simulations.to_string(simulation), particle_type_graph_names[j])),
                              cv2.VideoWriter_fourcc(*'DIVX'),
                              8,
                              frameSize)
    
        for filename in glob.glob(os.path.join(image_folder, "{} {} *.png".format(Simulations.to_string(simulation), particle_type_graph_names[j]))):
            img = cv2.imread(filename)
            out.write(img)
    
        out.release()
