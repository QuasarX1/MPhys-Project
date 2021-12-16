import numpy as np
from matplotlib import pyplot as plt

from DataAccess import constants, load_catalouge_field, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem

relitive_data_root = ".\\gm_for_mphys"

# GroupNumber (Subhalo) = index in total FOF table
# SubGroupNumber (Subhalo) = Subfind sub-group
# FirstSubhaloID (FOF) = index of first subhalo in Subhalo
# NumOfSubhalos (FOF) = number of values inc. first subhalo to read sequentially from Subhalo

"""
At each snapshot, working backwards, find the 100 most bound DARK MATTER particles in the identified object.
In the previous snapshot, identify the object containing the greatest number of most bound DARK MATTER particles from the previous itteration.
Save GroupNumber and SubGroupNumber.
Repeat.
"""

def get_mout_bound_particles(tag, simulation, halo, subhalo):
    n = 100
    group_particles = particles_by_group(tag, ParticleType.dark_matter, simulation, SimulationModels.RECAL, relitive_data_root)






centeral_mass_group = [[], [], []]
centeral_mass_sub_group = [[], [], []]
physical_centeral_mass_values = [[], [], []]
physical_centeral_mass_200_values = [[], [], []]
for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
    for j, tag in enumerate(constants.tags):
        if j == 0:
            continue

        try:
            object_masses = load_catalouge_field("GroupMass", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            object_masses_200 = load_catalouge_field("Group_M_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            #object_subhalo_IDs = load_catalouge_field("FirstSubhaloID", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            #group_numbers = load_catalouge_field("GroupNumber", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            #subgroup_numbers = load_catalouge_field("SubGroupNumber", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)

            field_selection_filter = object_masses == object_masses.max()
            #field_selection_filter = object_masses_200 == object_masses_200.max()
            physical_centeral_mass_values[i].append(object_masses[field_selection_filter][0])
            physical_centeral_mass_200_values[i].append(object_masses_200[field_selection_filter][0])
            
            #field_selection_filter = subgroup_numbers == object_subhalo_IDs[field_selection_filter][0]
            #centeral_mass_sub_group[i].append(subgroup_numbers[field_selection_filter][0])
            #centeral_mass_group[i].append(group_numbers[field_selection_filter][0])
        except LookupError:
            #centeral_mass_group[i].append(None)
            #centeral_mass_sub_group[i].append(None)
            physical_centeral_mass_values[i].append(None)

centeral_mass_group = np.array(centeral_mass_group)
centeral_mass_sub_group = np.array(centeral_mass_sub_group)
physical_centeral_mass_values = np.array(physical_centeral_mass_values)
physical_centeral_mass_200_values = np.array(physical_centeral_mass_200_values)



#for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
#    for j, tag in enumerate(constants.tags):
#        if j == 0:
#            continue
#
#        particles = particles_by_group(tag, ParticleType.gas, simulation, SimulationModels.RECAL, relitive_data_root)



snapshot_numbers = list(range(1, len(constants.tags)))
plt.plot(snapshot_numbers, physical_centeral_mass_200_values[0], label = "Early")
plt.plot(snapshot_numbers, physical_centeral_mass_200_values[1], label = "Organic")
plt.plot(snapshot_numbers, physical_centeral_mass_200_values[2], label = "Late")
plt.xlabel("Snapshot Number")
plt.ylabel("$M_{200}$")
plt.title("Group_M_Crit200 for each Snapshot")
plt.legend()
plt.show()

snapshot_numbers = list(range(1, len(constants.tags)))
plt.plot(snapshot_numbers, physical_centeral_mass_values[0], label = "Early")
plt.plot(snapshot_numbers, physical_centeral_mass_values[1], label = "Organic")
plt.plot(snapshot_numbers, physical_centeral_mass_values[2], label = "Late")
plt.xlabel("Snapshot Number")
plt.ylabel("Group Mass")
plt.title("GroupMass for each Snapshot")
plt.legend()
plt.show()
