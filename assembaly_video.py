import numpy as np
from matplotlib import pyplot as plt
import os
import glob
import cv2

from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, Simulations, SimulationModels, ParticleType, UnitSystem

from assembily_history import physical_centeral_mass_positions

model = SimulationModels.RECAL

relitive_data_root = ".\\gm_for_mphys"

image_folder = ".\\image_dump"
if not os.path.exists(image_folder):
    os.mkdir(image_folder)



particle_types = (ParticleType.gas, ParticleType.star, ParticleType.dark_matter)
particle_type_graph_names = ("Gas", "Star" ,"Dark Matter")

simulations = (Simulations.Early, Simulations.Organic, Simulations.Late)
for i in range(len(simulations)):
    for tag in constants.tags:
        try:
            object_masses = load_catalouge_field("GroupMass", "FOF", tag, simulations[i], SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            potential_centre = load_catalouge_field("GroupCentreOfPotential", "FOF", tag, simulations[i], SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)[object_masses == object_masses.max()][0]
            r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulations[i], SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)[object_masses == object_masses.max()][0]
            r_200 = r_200 / 8#TODO: why is the R_200 too large???
        except LookupError:
            continue# No objects in this snapshot - probably the first snapshot so just start the next one

        snapshot = ParticleReadConversion_EagleSnapshot(tag, simulations[i], model, relitive_data_root)

        for j in range(len(particle_types)):
            try:
                #particle_locations = snapshot.particle_read_sphere(particle_types[j], "Coordinates", potential_centre, r_200, UnitSystem.physical, UnitSystem.physical)

                galaxy_picture_edge_offsets = [0.15, 0.025, 0.025]
                lower_physical_limits = potential_centre + np.full((3, ), -galaxy_picture_edge_offsets[i])
                upper_physical_limits = potential_centre + np.full((3, ), galaxy_picture_edge_offsets[i])
                particle_locations = snapshot.centre_particles(snapshot.particle_read(particle_types[j], "Coordinates", lower_physical_limits, upper_physical_limits, UnitSystem.physical, UnitSystem.physical), potential_centre)
            except KeyError:
                continue# No particles of this type avalible

            fig = plt.figure(figsize = (12, 12))

            ax = plt.subplot(2, 2, 1)
            ax.set_xlabel("X (pMpc)")
            ax.set_ylabel("Y (pMpc)")
            ax.scatter(particle_locations[:, 0], particle_locations[:, 1], s = 0.01)

            ax = plt.subplot(2, 2, 2)
            ax.set_xlabel("Z (pMpc)")
            ax.set_ylabel("Y (pMpc)")
            ax.scatter(particle_locations[:, 2], particle_locations[:, 1], s = 0.01)

            ax = plt.subplot(2, 2, 3)
            ax.set_xlabel("X (pMpc)")
            ax.set_ylabel("Z (pMpc)")
            ax.scatter(particle_locations[:, 0], particle_locations[:, 2], s = 0.01)
    
            fig.suptitle("{} Particle Projection of the Assembled Galaxy in the {} Regime".format(particle_type_graph_names[j], Simulations.to_string(simulations[i]).title()))
            plt.savefig(os.path.join(image_folder, "{} {} {}.png".format(Simulations.to_string(simulations[i]).title(), particle_type_graph_names[j], tag)))

    # Create Video
    frameSize = (1200, 1200)
    for j in range(len(particle_types)):
        out = cv2.VideoWriter(os.path.join(image_folder, "{} {}.avi".format(Simulations.to_string(simulations[i]), particle_type_graph_names[j])),
                              cv2.VideoWriter_fourcc(*'DIVX'),
                              8,
                              frameSize)
    
        for filename in glob.glob(os.path.join(image_folder, "{} {} *.png".format(Simulations.to_string(simulations[i]), particle_type_graph_names[j]))):
            img = cv2.imread(filename)
            out.write(img)
    
        out.release()
