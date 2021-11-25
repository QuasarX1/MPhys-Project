import os
import numpy as np
from matplotlib import pyplot as plt

import h5py as h5
import pyread_eagle

from DataAccess import ParticleReadConversion_EagleSnapshot, Simulations, SimulationModels, ParticleType, UnitSystem

def combine_limits(lower, upper):
    bound_length = len(lower)
    if bound_length != len(upper):
        raise ValueError("Arrays of bounds were not of the same length.")
    return np.append(np.array(lower).reshape((bound_length, 1)), np.array(upper).reshape((bound_length, 1)), axis = 1).reshape(2 * bound_length)

sim = Simulations.Organic
model = SimulationModels.RECAL
tag = "028_z000p000"

relitive_data_root = ".\\gm_for_mphys"

snapshot = ParticleReadConversion_EagleSnapshot(tag, sim, model, relitive_data_root)

print("Box size in GADGET units = {}".format(snapshot.boxsize))
print("Box size in physical units = {}".format(snapshot.convert_distance_values(snapshot.boxsize, UnitSystem.h_less_comoving_GADGET, UnitSystem.physical)))
print()

galaxy_centre = np.array([28.326, 2.614, 21.519])
galaxy_bounds_lower_relitive = np.array([-0.025, -0.025, -0.025])
galaxy_bounds_upper_relitive = np.array([0.025, 0.025, 0.025])

galaxy_centre_physical_units = snapshot.convert_distance_values(galaxy_centre, UnitSystem.h_less_comoving_GADGET, UnitSystem.physical)
lower_physical_limits = galaxy_bounds_lower_relitive + galaxy_centre_physical_units
upper_physical_limits = galaxy_bounds_upper_relitive + galaxy_centre_physical_units

lower_gadget_limits = snapshot.convert_distance_values(lower_physical_limits, UnitSystem.physical, UnitSystem.h_less_comoving_GADGET)
upper_gadget_limits = snapshot.convert_distance_values(upper_physical_limits, UnitSystem.physical, UnitSystem.h_less_comoving_GADGET)


#limits[0.0, snapshot.boxsize, 0.0, snapshot.boxsize, 0.0, snapshot.boxsize]
limits = combine_limits(lower_gadget_limits, upper_gadget_limits)

print("Region bounds in GADGET units = {}".format(limits))
print("Region bounds in relitive physical units = {}".format(combine_limits(galaxy_bounds_lower_relitive, galaxy_bounds_upper_relitive)))
print()

snapshot.select_region(*limits)






#gas_particle_locations = snapshot.particle_read(ParticleType.gas, "Coordinates", UnitSystem.physical)
gas_particle_locations = snapshot.convert_distance_values(snapshot.particle_read(ParticleType.gas, "Coordinates"), UnitSystem.h_less_comoving_GADGET, UnitSystem.physical)
#print(gas_particle_locations)
#print()

gas_particle_locations_centre_adjusted = gas_particle_locations - galaxy_centre_physical_units
#print(gas_particle_locations_centre_adjusted)
#print()

print("X Data Extremes:")
print(min(gas_particle_locations_centre_adjusted[:, 0]))
print(max(gas_particle_locations_centre_adjusted[:, 0]))
print("Y Data Extremes:")
print(min(gas_particle_locations_centre_adjusted[:, 1]))
print(max(gas_particle_locations_centre_adjusted[:, 1]))
print("Z Data Extremes:")
print(min(gas_particle_locations_centre_adjusted[:, 2]))
print(max(gas_particle_locations_centre_adjusted[:, 2]))
print()

#plt.scatter(gas_particle_locations_centre_adjusted[:, 0], gas_particle_locations_centre_adjusted[:, 1])
#plt.show()
ax = plt.subplot(2, 2, 1)
ax.scatter(gas_particle_locations_centre_adjusted[:, 0], gas_particle_locations_centre_adjusted[:, 1], s = 0.01)

#plt.scatter(gas_particle_locations_centre_adjusted[:, 2], gas_particle_locations_centre_adjusted[:, 1])
#plt.show()
ax = plt.subplot(2, 2, 2)
ax.scatter(gas_particle_locations_centre_adjusted[:, 2], gas_particle_locations_centre_adjusted[:, 1], s = 0.01)

#plt.scatter(gas_particle_locations_centre_adjusted[:, 0], gas_particle_locations_centre_adjusted[:, 2])
#plt.show()
ax = plt.subplot(2, 2, 3)
ax.scatter(gas_particle_locations_centre_adjusted[:, 0], gas_particle_locations_centre_adjusted[:, 2], s = 0.01)

plt.show()

# x = (-0.025, 0.025), y = (-0.025, 0.025), z = (-0.025, 0.025)
