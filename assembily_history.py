import numpy as np

from DataAccess import constants, load_catalouge_field, Simulations, SimulationModels, ParticleType, UnitSystem

relitive_data_root = ".\\gm_for_mphys"

physical_centeral_mass_positions = [[], [], []]
for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
    for tag in constants.tags:
        try:
            object_masses = load_catalouge_field("GroupMass", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            potential_centre = load_catalouge_field("GroupCentreOfPotential", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.physical)
            physical_centeral_mass_positions[i].append(potential_centre[object_masses == object_masses.max()][0])
        except LookupError:
            physical_centeral_mass_positions[i].append(None)

physical_centeral_mass_positions = np.array(physical_centeral_mass_positions)

#print(physical_centeral_mass_positions[0, -1] - physical_centeral_mass_positions[1, -1])
#print(physical_centeral_mass_positions[1, -1])
#print(physical_centeral_mass_positions[2, -1] - physical_centeral_mass_positions[1, -1])
