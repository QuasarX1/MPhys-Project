import numpy as np
from matplotlib import pyplot as plt

from DataAccess import ParticleReadConversion_EagleSnapshot, Simulations, SimulationModels, ParticleType, UnitSystem, constants

from assembily_history import physical_centeral_mass_positions

sim = Simulations.Organic
model = SimulationModels.RECAL
tag = "028_z000p000"

relitive_data_root = ".\\gm_for_mphys"

#snapshot = ParticleReadConversion_EagleSnapshot(tag, sim, model, relitive_data_root)

#physical_box_size = snapshot.convert_distance_values(snapshot.boxsize, UnitSystem.h_less_comoving_GADGET, UnitSystem.physical)
#print("Box size in GADGET units = {}".format(snapshot.boxsize))
#print("Box size in physical units = {}".format(physical_box_size))
#print()

#organic_galaxy_centre = np.array([28.326, 2.614, 21.519])
#organic_galaxy_centre_physical_units = snapshot.convert_distance_values(organic_galaxy_centre, UnitSystem.h_less_comoving_GADGET, UnitSystem.physical)

#galaxy_centres = [organic_galaxy_centre_physical_units + np.array([0.42436218, 0.63016963, 0.34677315]),
#                  organic_galaxy_centre_physical_units,
#                  organic_galaxy_centre_physical_units + np.array([-0.839489, -0.5783951, -0.05172157])]

#galaxy_centres = physical_centeral_mass_positions[:, -1]
galaxy_centres = [physical_centeral_mass_positions[sim][constants.tags[-1]] for sim in Simulations]

galaxy_picture_edge_offsets = [0.15, 0.025, 0.025]



"""
snapshot = ParticleReadConversion_EagleSnapshot(tag, Simulations.Early, model, relitive_data_root)

offset = galaxy_picture_edge_offsets[0]
galaxy_centre_physical_units = galaxy_centres[0]

galaxy_bounds_lower_relitive = np.full((3, ), -offset)
galaxy_bounds_upper_relitive = np.full((3, ), offset)

lower_physical_limits = galaxy_bounds_lower_relitive + galaxy_centre_physical_units
upper_physical_limits = galaxy_bounds_upper_relitive + galaxy_centre_physical_units

#particle_locations = snapshot.particle_read(ParticleType.star, "Coordinates", unit_system = UnitSystem.physical)
particle_locations = snapshot.particle_read(ParticleType.star, "Coordinates", lower_physical_limits, upper_physical_limits, UnitSystem.physical, UnitSystem.physical)

particle_locations_centre_adjusted = particle_locations - galaxy_centre_physical_units

particle_locations_box_adjusted = ((particle_locations_centre_adjusted + (physical_box_size / 2)) % physical_box_size) - (physical_box_size / 2)

fig = plt.figure(figsize = (12, 12))

ax = plt.subplot(2, 2, 1)
ax.set_xlabel("X (pMpc)")
ax.set_ylabel("Y (pMpc)")
ax.scatter(particle_locations_box_adjusted[:, 0], particle_locations_box_adjusted[:, 1], s = 0.01)

ax = plt.subplot(2, 2, 2)
ax.set_xlabel("Z (pMpc)")
ax.set_ylabel("Y (pMpc)")
ax.scatter(particle_locations_box_adjusted[:, 2], particle_locations_box_adjusted[:, 1], s = 0.01)

ax = plt.subplot(2, 2, 3)
ax.set_xlabel("X (pMpc)")
ax.set_ylabel("Z (pMpc)")
ax.scatter(particle_locations_box_adjusted[:, 0], particle_locations_box_adjusted[:, 2], s = 0.01)

plt.show()
"""



simulations = (Simulations.Early, Simulations.Organic, Simulations.Late)
for i in range(len(simulations)):
    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulations[i], model, relitive_data_root)

    offset = galaxy_picture_edge_offsets[i]
    galaxy_centre_physical_units = galaxy_centres[i]

    galaxy_bounds_lower_relitive = np.full((3, ), -offset)
    galaxy_bounds_upper_relitive = np.full((3, ), offset)
    
    lower_physical_limits = galaxy_bounds_lower_relitive + galaxy_centre_physical_units
    upper_physical_limits = galaxy_bounds_upper_relitive + galaxy_centre_physical_units

    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.dark_matter)
    particle_type_graph_names = ("Gas", "Star" ,"Dark Matter")
    for j in range(len(particle_types)):
        particle_locations = snapshot.particle_read(particle_types[j], "Coordinates", lower_physical_limits, upper_physical_limits, UnitSystem.physical, UnitSystem.physical)

        #particle_locations_centre_adjusted = particle_locations - galaxy_centre_physical_units

        #particle_locations_box_adjusted = ((particle_locations_centre_adjusted + (physical_box_size / 2)) % physical_box_size) - (physical_box_size / 2)

        particle_locations_box_adjusted = snapshot.centre_particles(particle_locations, galaxy_centre_physical_units)

        fig = plt.figure(figsize = (12, 12))

        ax = plt.subplot(2, 2, 1)
        ax.set_xlabel("X (pMpc)")
        ax.set_ylabel("Y (pMpc)")
        ax.scatter(particle_locations_box_adjusted[:, 0], particle_locations_box_adjusted[:, 1], s = 0.01)

        ax = plt.subplot(2, 2, 2)
        ax.set_xlabel("Z (pMpc)")
        ax.set_ylabel("Y (pMpc)")
        ax.scatter(particle_locations_box_adjusted[:, 2], particle_locations_box_adjusted[:, 1], s = 0.01)

        ax = plt.subplot(2, 2, 3)
        ax.set_xlabel("X (pMpc)")
        ax.set_ylabel("Z (pMpc)")
        ax.scatter(particle_locations_box_adjusted[:, 0], particle_locations_box_adjusted[:, 2], s = 0.01)
    
        fig.suptitle("{} Particle Projection of the Assembled Galaxy in the {} Regime".format(particle_type_graph_names[j], Simulations.to_string(simulations[i]).title()))
        plt.show()
"""
"""

# x = 0.42 (-0.15, 0.15), y = 0.63 (-0.15, 0.15), z = 0.34 (-0.15, 0.15)
# x = (-0.025, 0.025), y = (-0.025, 0.025), z = (-0.025, 0.025)
