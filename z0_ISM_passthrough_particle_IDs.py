from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, UnitSystem
from assembily_history import physical_centeral_mass_positions, assembily_history
import os
import numpy as np

relitive_data_root = ".\\gm_for_mphys"

def find_particle_data(sim: Simulations, r_200_fraction: float = 1.0):
    tag = constants.tags[-1]
    particleType = ParticleType.gas

    snapshot = ParticleReadConversion_EagleSnapshot(tag, sim, SimulationModels.RECAL, relitive_data_root)

    if r_200_fraction is not None:
        r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][tag]["halo"] - 1]
        selection_radius = r_200 * r_200_fraction

        galaxy_centre = physical_centeral_mass_positions[sim][tag]# cgs
        
        IDs = snapshot.particle_read_sphere(particleType, "ParticleIDs", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        star_forming_history_metric = snapshot.particle_read_sphere(particleType, "OnEquationOfState", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        max_temp = snapshot.particle_read_sphere(particleType, "MaximumTemperature", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        expansion_factor_at_max_temp = snapshot.particle_read_sphere(particleType, "AExpMaximumTemperature", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)

    else:
        IDs = snapshot.particle_read(particleType, "ParticleIDs", unit_system = UnitSystem.cgs)
        star_forming_history_metric = snapshot.particle_read(particleType, "OnEquationOfState", unit_system = UnitSystem.cgs)
        max_temp = snapshot.particle_read(particleType, "MaximumTemperature", unit_system = UnitSystem.cgs)
        expansion_factor_at_max_temp = snapshot.particle_read(particleType, "AExpMaximumTemperature", unit_system = UnitSystem.cgs)

    selectionFilter = star_forming_history_metric < 0

    IDs = IDs[selectionFilter]
    star_forming_history_metric = star_forming_history_metric[selectionFilter]
    lastStarFormingAt = star_forming_history_metric * -1# Expansion factor when particle was last star forming
    max_temp = max_temp[selectionFilter]
    expansion_factor_at_max_temp = expansion_factor_at_max_temp[selectionFilter]

    with open("{}_final_ISM_{}_particle_IDs.txt".format(Simulations.to_string(sim), "R_200_times_{}".format(r_200_fraction) if r_200_fraction is not None else "whole_box"), "w") as f:
        f.write("--|| PARTICLE TYPE ||--:{}\n".format(ParticleType.to_string(particleType)))
        for i in range(len(IDs)):
            f.write("{}:{}:{}:{}\n".format(IDs[i], lastStarFormingAt[i], max_temp[i], expansion_factor_at_max_temp[i]))
    
    

def check_exists_patricle_data(sim: Simulations, r_200_fraction: float = 1.0):
    return os.path.exists("{}_final_ISM_{}_particle_IDs.txt".format(Simulations.to_string(sim), "R_200_times_{}".format(r_200_fraction) if r_200_fraction is not None else "whole_box"))



def load_particle_data(sim: Simulations, r_200_fraction: float = 1.0):
    if not check_exists_patricle_data(sim, r_200_fraction):
        find_particle_data(sim, r_200_fraction)
        
    particle_IDs_by_type = {}

    with open("{}_final_ISM_{}_particle_IDs.txt".format(Simulations.to_string(sim), "R_200_times_{}".format(r_200_fraction) if r_200_fraction is not None else "whole_box"), "r") as f:
        lines = f.readlines()
        IDs = []
        expansion_factors_last_star_forming = []
        max_temps = []
        expansion_factors_max_temp = []
        for line in lines:
            line = line[:-1]# Remove the newline character
            if line.split(":", 1)[0] == "--|| PARTICLE TYPE ||--":
                pass
            else:
                ID, a_last_star_forming, max_temp, a_max_temp = line.split(":")
                IDs.append(ID)
                expansion_factors_last_star_forming.append(a_last_star_forming)
                max_temps.append(max_temp)
                expansion_factors_max_temp.append(a_max_temp)


    return np.array(IDs, dtype = float), np.array(expansion_factors_last_star_forming, dtype = float), np.array(max_temps, dtype = float), np.array(expansion_factors_max_temp, dtype = float)
