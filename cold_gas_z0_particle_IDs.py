from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, UnitSystem
from assembily_history import physical_centeral_mass_positions, assembily_history
import os

relitive_data_root = ".\\gm_for_mphys"

def find_particle_IDs(sim: Simulations, r_200_fraction: float = 1.0):
    particle_IDs_by_type = {}

    tag = constants.tags[-1]
    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][tag]["halo"] - 1]
    
    selection_radius = r_200 * r_200_fraction

    galaxy_centre = physical_centeral_mass_positions[sim][tag]# cgs
    
    snapshot = ParticleReadConversion_EagleSnapshot(tag, sim, SimulationModels.RECAL, relitive_data_root)
    for particle_type in all_particle_types:
        particles = snapshot.particle_read_sphere(particle_type, "ParticleIDs", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        if particle_type == ParticleType.gas:
            formation_rates = snapshot.particle_read_sphere(particle_type, "StarFormationRate", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
            particles = particles[formation_rates == 0]# Cold gas is not star forming
        particle_IDs_by_type[particle_type] = particles
        
    with open("{}_final_R_200_times_{}_particle_IDs_cold_gas_only.txt".format(Simulations.to_string(sim), r_200_fraction), "w") as f:
        for particle_type in all_particle_types:
            f.write("--|| PARTICLE TYPE ||--:{}\n".format(ParticleType.to_string(particle_type)))
            for value in particle_IDs_by_type[particle_type]:
                f.write("{}\n".format(value))
    
    

def check_exists_patricle_IDs(sim: Simulations, r_200_fraction: float = 1.0):
    return os.path.exists("{}_final_R_200_times_{}_particle_IDs_cold_gas_only.txt".format(Simulations.to_string(sim), r_200_fraction))



def load_particle_IDs(sim: Simulations, r_200_fraction: float = 1.0):
    if not check_exists_patricle_IDs(sim, r_200_fraction):
        find_particle_IDs(sim, r_200_fraction)
        
    particle_IDs_by_type = {}

    with open("{}_final_R_200_times_{}_particle_IDs_cold_gas_only.txt".format(Simulations.to_string(sim), r_200_fraction), "r") as f:
        lines = f.readlines()
        current_particle_type = None
        for line in lines:
            line = line[:-1]# Remove the newline character
            if line.split(":", 1)[0] == "--|| PARTICLE TYPE ||--":
                current_particle_type = ParticleType.from_string(line.split(":", 1)[1])
                particle_IDs_by_type[current_particle_type] = []
            else:
                particle_IDs_by_type[current_particle_type].append(line)


    return particle_IDs_by_type



if __name__ == "__main__":
    if not check_exists_patricle_IDs(Simulations.Early, 1.0):
        find_particle_IDs(Simulations.Early, 1.0)
    
    if not check_exists_patricle_IDs(Simulations.Organic, 1.0):
        find_particle_IDs(Simulations.Organic, 1.0)
    
    if not check_exists_patricle_IDs(Simulations.Late, 1.0):
        find_particle_IDs(Simulations.Late, 1.0)