import numpy as np
import h5py as h5
import pyread_eagle
import os
from matplotlib import pyplot as plt

from DataAccess import ParticleReadConversion_EagleSnapshot, Simulations, SimulationModels, ParticleType, UnitSystem

sim = "Organic"
model = "RECAL"
tag = "019_z001p004"

relitive_data_root = "..\\..\\gm_for_mphys"



#snapfile_template = os.path.join(relitive_data_root, model, sim, f"snapshot_{tag}", f"snap_{tag}.{{}}.hdf5")
#snapfile = snapfile_template.format(0)
#
#snapshot = pyread_eagle.EagleSnapshot(snapfile)
#
#with h5.File(snapfile,'r') as datafile:
#    h = datafile['Header'].attrs['HubbleParam']
#    a = datafile['Header'].attrs['ExpansionFactor']
#
#side_length = 100.0 * h / a# Side length 1 pMpc



alternate_snapshot_object = ParticleReadConversion_EagleSnapshot("019_z001p004", Simulations.Organic, SimulationModels.RECAL, relitive_data_root)

print(alternate_snapshot_object.header)
side_length = 100.0 * alternate_snapshot_object.hubble_paramiter / alternate_snapshot_object.expansion_factor# Side length 1 pMpc

alternate_snapshot_object.select_region(0.0, side_length, 0.0, side_length, 0.0, side_length)
gas_mass_alt = alternate_snapshot_object.particle_read(ParticleType.gas, "Mass", UnitSystem.physical)
print("Alternate Gas Mass -", len(gas_mass_alt), "items:", gas_mass_alt)

gas_locations = alternate_snapshot_object.particle_read(ParticleType.gas, "Coordinates", UnitSystem.physical)
plt.scatter(gas_locations[:, 0], gas_locations[:, 1])
plt.show()
plt.scatter(gas_locations[:, 2], gas_locations[:, 1])
plt.show()
plt.scatter(gas_locations[:, 0], gas_locations[:, 2])
plt.show()
