import numpy as np

from DataAccess import constants, load_catalouge_field, Simulations, SimulationModels, ParticleType, UnitSystem

relitive_data_root = ".\\gm_for_mphys"

assembily_history = {
    Simulations.Early:{
        "000_z020p000": {"halo":None, "subhalo":None},
        "001_z015p132": {"halo":None, "subhalo":None},
        "002_z009p993": {"halo":7,    "subhalo":0   },
        "003_z008p988": {"halo":5,    "subhalo":0   },
        "004_z008p075": {"halo":4,    "subhalo":0   },
        "005_z007p050": {"halo":1,    "subhalo":2   },
        "006_z005p971": {"halo":2,    "subhalo":1   },
        "007_z005p487": {"halo":1,    "subhalo":0   },
        "008_z005p037": {"halo":2,    "subhalo":0   },
        "009_z004p485": {"halo":4,    "subhalo":0   },
        "010_z003p984": {"halo":1,    "subhalo":1   },
        "011_z003p528": {"halo":1,    "subhalo":0   },
        "012_z003p017": {"halo":1,    "subhalo":0   },
        "013_z002p478": {"halo":1,    "subhalo":0   },
        "014_z002p237": {"halo":1,    "subhalo":0   },
        "015_z002p012": {"halo":1,    "subhalo":0   },
        "016_z001p737": {"halo":None, "subhalo":None},
        "017_z001p487": {"halo":1,    "subhalo":0   },
        "018_z001p259": {"halo":1,    "subhalo":0   },
        "019_z001p004": {"halo":1,    "subhalo":0   },
        "020_z000p865": {"halo":1,    "subhalo":0   },
        "021_z000p736": {"halo":1,    "subhalo":0   },
        "022_z000p615": {"halo":1,    "subhalo":0   },
        "023_z000p503": {"halo":1,    "subhalo":0   },
        "024_z000p366": {"halo":1,    "subhalo":0   },
        "025_z000p271": {"halo":1,    "subhalo":0   },
        "026_z000p183": {"halo":1,    "subhalo":0   },
        "027_z000p101": {"halo":1,    "subhalo":0   },
        "028_z000p000": {"halo":1,    "subhalo":0   }
    },

    Simulations.Organic:{
        "000_z020p000": {"halo":None, "subhalo":None},
        "001_z015p132": {"halo":None, "subhalo":None},
        "002_z009p993": {"halo":6,    "subhalo":0   },
        "003_z008p988": {"halo":14,   "subhalo":0   },
        "004_z008p075": {"halo":7,    "subhalo":0   },
        "005_z007p050": {"halo":8,    "subhalo":0   },
        "006_z005p971": {"halo":11,   "subhalo":0   },
        "007_z005p487": {"halo":9,    "subhalo":0   },
        "008_z005p037": {"halo":3,    "subhalo":1   },
        "009_z004p485": {"halo":3,    "subhalo":0   },
        "010_z003p984": {"halo":4,    "subhalo":0   },
        "011_z003p528": {"halo":5,    "subhalo":0   },
        "012_z003p017": {"halo":5,    "subhalo":0   },
        "013_z002p478": {"halo":1,    "subhalo":0   },
        "014_z002p237": {"halo":1,    "subhalo":0   },
        "015_z002p012": {"halo":1,    "subhalo":0   },
        "016_z001p737": {"halo":1,    "subhalo":0   },
        "017_z001p487": {"halo":None, "subhalo":None},
        "018_z001p259": {"halo":1,    "subhalo":0   },
        "019_z001p004": {"halo":1,    "subhalo":0   },
        "020_z000p865": {"halo":1,    "subhalo":0   },
        "021_z000p736": {"halo":1,    "subhalo":0   },
        "022_z000p615": {"halo":1,    "subhalo":0   },
        "023_z000p503": {"halo":1,    "subhalo":0   },
        "024_z000p366": {"halo":1,    "subhalo":0   },
        "025_z000p271": {"halo":1,    "subhalo":0   },
        "026_z000p183": {"halo":1,    "subhalo":0   },
        "027_z000p101": {"halo":1,    "subhalo":0   },
        "028_z000p000": {"halo":1,    "subhalo":0   }
    },

    Simulations.Late:{
        "000_z020p000": {"halo":None, "subhalo":None},
        "001_z015p132": {"halo":3,    "subhalo":0   },
        "002_z009p993": {"halo":1,    "subhalo":0   },
        "003_z008p988": {"halo":1,    "subhalo":0   },
        "004_z008p075": {"halo":1,    "subhalo":0   },
        "005_z007p050": {"halo":1,    "subhalo":0   },
        "006_z005p971": {"halo":1,    "subhalo":0   },
        "007_z005p487": {"halo":1,    "subhalo":0   },
        "008_z005p037": {"halo":1,    "subhalo":0   },
        "009_z004p485": {"halo":1,    "subhalo":0   },
        "010_z003p984": {"halo":1,    "subhalo":0   },
        "011_z003p528": {"halo":1,    "subhalo":0   },
        "012_z003p017": {"halo":1,    "subhalo":0   },
        "013_z002p478": {"halo":1,    "subhalo":0   },
        "014_z002p237": {"halo":2,    "subhalo":0   },
        "015_z002p012": {"halo":1,    "subhalo":0   },
        "016_z001p737": {"halo":1,    "subhalo":0   },
        "017_z001p487": {"halo":1,    "subhalo":0   },
        "018_z001p259": {"halo":1,    "subhalo":0   },
        "019_z001p004": {"halo":1,    "subhalo":0   },
        "020_z000p865": {"halo":1,    "subhalo":0   },
        "021_z000p736": {"halo":1,    "subhalo":0   },
        "022_z000p615": {"halo":1,    "subhalo":0   },
        "023_z000p503": {"halo":1,    "subhalo":0   },
        "024_z000p366": {"halo":1,    "subhalo":0   },
        "025_z000p271": {"halo":1,    "subhalo":0   },
        "026_z000p183": {"halo":1,    "subhalo":0   },
        "027_z000p101": {"halo":1,    "subhalo":0   },
        "028_z000p000": {"halo":1,    "subhalo":0   }
    }
}

physical_centeral_mass_positions = {}
for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
#for i, simulation in enumerate((Simulations.Organic, )):#TODO: enable to only use organic
    physical_centeral_mass_positions[simulation] = {}
    for tag in constants.tags:
#    for tag in (constants.tags[-1], ):#TODO: enable to only use z=0
        try:
            halo = assembily_history[Simulations.Early][tag]["halo"]
            subhalo = assembily_history[Simulations.Early][tag]["subhalo"]
            
            if halo is None or subhalo is None:
                raise LookupError("No history avalible.")

            #potential_centre = load_catalouge_field("CentreOfPotential", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            potential_centre = load_catalouge_field("CentreOfPotential", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)
            group_number = load_catalouge_field("GroupNumber", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root)
            sub_group_number = load_catalouge_field("SubGroupNumber", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root)
            physical_centeral_mass_positions[simulation][tag] = potential_centre[(group_number == halo) & (sub_group_number == subhalo)][0]
        except LookupError:
            physical_centeral_mass_positions[simulation][tag] = None



"""
import numpy as np
from matplotlib import pyplot as plt

from DataAccess import constants, load_catalouge_field, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem

relitive_data_root = ".\\gm_for_mphys"

# GroupNumber (Subhalo) = index in total FOF table
# SubGroupNumber (Subhalo) = Subfind sub-group
# FirstSubhaloID (FOF) = index of first subhalo in Subhalo
# NumOfSubhalos (FOF) = number of values inc. first subhalo to read sequentially from Subhalo

'''
At each snapshot, working backwards, find the 100 most bound DARK MATTER particles in the identified object.
In the previous snapshot, identify the object containing the greatest number of most bound DARK MATTER particles from the previous itteration.
Save GroupNumber and SubGroupNumber.
Repeat.
'''

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
