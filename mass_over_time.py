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

def get_most_bound_particles(tag, simulation, halo, subhalo):
    n = 100
    fields = ["ParticleIDs", "GroupNumber", "SubGroupNumber", "ParticleBindingEnergy"]
    print(f"Reading data for tag {tag}:")
    group_particles = particles_by_group(tag, ParticleType.dark_matter, simulation, SimulationModels.RECAL, relitive_data_root, fields, verbose = True)#TODO: remove verbose
    print()
    group_particles = group_particles[:, (group_particles[1, :] == halo) & (group_particles[2, :] == subhalo)]
    
    most_bound_IDs = []
    
    most_bound_IDs.append(group_particles[0, np.where(group_particles[3, :] == group_particles[3, :].min())[0]][0])
    for _ in range(1, n):
        most_bound_IDs.append(group_particles[:, True ^ np.in1d(group_particles[0, :], most_bound_IDs)][0, np.where(group_particles[3, :] == group_particles[3, :].min())[0]][0])

    return most_bound_IDs

def get_greatest_particle_membership(tag, simulation, particle_IDs):
    fields = ["ParticleIDs", "GroupNumber", "SubGroupNumber"]
    all_particles = particles_by_group(tag, ParticleType.dark_matter, simulation, SimulationModels.RECAL, relitive_data_root, fields)
    selected_particles = all_particles[:, np.in1d(all_particles[0, :], particle_IDs)]

    halo_values, halo_counts = np.unique(selected_particles[1, :], return_counts=True)
    if len(halo_counts) > 0:
        most_common_halo = halo_values[halo_counts.argmax()]

        sub_halo_values, sub_halo_counts = np.unique(selected_particles[:, selected_particles[1, :] == most_common_halo][2, :], return_counts=True)
        most_common_sub_halo = sub_halo_values[sub_halo_counts.argmax()] if len(sub_halo_counts) > 0 else None

    else:# None of the requested particles were part of a registered halo or subhalo
        most_common_halo = None
        most_common_sub_halo = None

    return most_common_halo, most_common_sub_halo



halo_numbers = [[], [], []]
sub_halo_numbers = [[], [], []]
reversed_tags = list(constants.tags)
reversed_tags.reverse()
#for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
for i, simulation in enumerate((Simulations.Organic, Simulations.Organic, Simulations.Organic)):
    max_mass_index = load_catalouge_field("Mass", "Subhalo", constants.tags[-1], simulation, SimulationModels.RECAL, relitive_data_root).argmax()
    current_halo = load_catalouge_field("GroupNumber", "Subhalo", constants.tags[-1], simulation, SimulationModels.RECAL, relitive_data_root)[max_mass_index]
    current_sub_halo = load_catalouge_field("SubGroupNumber", "Subhalo", constants.tags[-1], simulation, SimulationModels.RECAL, relitive_data_root)[max_mass_index]

    
    skipFile = False
    for j, tag in enumerate(reversed_tags):
        if not skipFile:
            halo_numbers[i].append(current_halo)
            sub_halo_numbers[i].append(current_sub_halo)
        else:
            halo_numbers[i].append(None)
            sub_halo_numbers[i].append(None)

        if j == len(reversed_tags) - 2:
            break

        if not skipFile:
            particle_IDs = get_most_bound_particles(tag, simulation, current_halo, current_sub_halo)

        skipFile = False

        try:
            current_halo, current_sub_halo = get_greatest_particle_membership(reversed_tags[j + 1], simulation, particle_IDs)
        except OSError:
            skipFile = True
        if current_sub_halo == None:
            skipFile = True

    halo_numbers[i].append(None)
    sub_halo_numbers[i].append(None)

    halo_numbers[i].reverse()
    sub_halo_numbers[i].reverse()

halo_numbers = np.array(halo_numbers)
sub_halo_numbers = np.array(sub_halo_numbers)

for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
    print(Simulations.to_string(simulation))
    for j in range(len(constants.tags)):
        print(constants.tags[j], halo_numbers[i][j], sub_halo_numbers[i][j])
    print()


"""
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
"""
