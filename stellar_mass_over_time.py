import numpy as np
from matplotlib import pyplot as plt

from DataAccess import constants, load_catalouge_field, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem
from assembily_history import assembily_history

relitive_data_root = ".\\gm_for_mphys"

mass_values = []
snapshot_numbers = []
for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
    mass_values.append(np.empty(shape = (len(constants.tags), )))
    snapshot_numbers.append(np.array([float("{}.{}".format(tag[5:8], tag[9:12])) for tag in constants.tags]))
    #snapshot_numbers.append(np.array(constants.tags))
    for j, tag in enumerate(constants.tags):
        try:
            halo = assembily_history[simulation][tag]["halo"]
            subhalo = assembily_history[simulation][tag]["subhalo"]
            if halo is None:
                raise LookupError("No history avalible.")

            first_subhalo_index = load_catalouge_field("FirstSubhaloID", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)[halo - 1]
            subhalo_index = first_subhalo_index + subhalo
            mass_values[i][j] = load_catalouge_field("ApertureMeasurements/Mass/030kpc", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)[subhalo_index][ParticleType.star.value]

        except LookupError:
            mass_values[i][j] = None
            
    array_filter = True ^ np.isnan(mass_values[i])
    mass_values[i] = mass_values[i][array_filter]
    snapshot_numbers[i] = snapshot_numbers[i][array_filter]




plt.plot(snapshot_numbers[0], mass_values[0], label = "Early")
plt.plot(snapshot_numbers[1], mass_values[1], label = "Organic")
plt.plot(snapshot_numbers[2], mass_values[2], label = "Late")
plt.loglog()
#plt.semilogx()
plt.xlim(plt.xlim()[1], plt.xlim()[0])
plt.xlabel("Redshift")
plt.ylabel("M*")
plt.title("Central Subhalo Stellar Mass for each Snapshot")
plt.legend()
plt.show()
