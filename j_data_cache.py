import numpy as np
from matplotlib import pyplot as plt
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, baryonic_particle_types, baryonic_low_mass_particle_types, UnitSystem
from assembily_history import physical_centeral_mass_positions, assembily_history
from z0_particle_IDs import load_particle_IDs
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs, X_Axis_Value, line_colours, line_styles
#from Physics import specific_angular_momentum, centre_of_mass
#from angular_momentum_units import specific_angular_momentum_units

from regional_specific_angular_momentum import set_particle_types, set_R_200_fraction, sphere_specific_angular_momentum, sphere_radius

relitive_data_root = ".\\gm_for_mphys"



if __name__ == "__main__":
    r_values_tags = {         ParticleType.gas : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                                   0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                             ParticleType.star : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                                   0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                      ParticleType.dark_matter : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                                   0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } }
                     }
    r_values = {         ParticleType.gas : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                              0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                        ParticleType.star : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                              0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                 ParticleType.dark_matter : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                              0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } }
                }

    j_values_tags = {         ParticleType.gas : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                                   0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                             ParticleType.star : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                                   0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                      ParticleType.dark_matter : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                                   0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } }
                     }
    j_values = {         ParticleType.gas : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                              0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                        ParticleType.star : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                              0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } },
                 ParticleType.dark_matter : { 1.0 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] },
                                              0.1 : { Simulations.Early : [], Simulations.Organic : [], Simulations.Late : [] } }
                }



    for particle_type in (ParticleType.gas, ParticleType.star, ParticleType.dark_matter):
        for r_200_fraction in (1.0, 0.1):
            for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
                for tag in constants.tags:
                    halo = assembily_history[sim][tag]["halo"]
                    subhalo = assembily_history[sim][tag]["subhalo"]
                    if halo is None:
                        continue
                    try:
                        r = sphere_radius(halo, subhalo, tag, sim, [particle_type], r_200_fraction)
                    except LookupError:
                        continue
                    r_values_tags[particle_type][r_200_fraction][sim].append(tag)
                    r_values[particle_type][r_200_fraction][sim].append(r)

    with open("r_values.txt", "w") as f:
        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
            f.write("--|| STYPE ||-- {}\n".format(Simulations.to_string(sim)))
            for r_200_fraction in (1.0, 0.1):
                f.write(f"--|| RATIO ||-- {r_200_fraction}\n")
                for particle_type in (ParticleType.gas, ParticleType.star, ParticleType.dark_matter):
                    f.write("--|| PTYPE ||-- {}\n".format(ParticleType.to_string(particle_type)))
                    for i in range(len(r_values_tags[particle_type][r_200_fraction][sim])):
                        f.write(f"{r_values_tags[particle_type][r_200_fraction][sim][i]}:{r_values[particle_type][r_200_fraction][sim][i]}\n")



    for particle_type in (ParticleType.gas, ParticleType.star, ParticleType.dark_matter):
        for r_200_fraction in (1.0, 0.1):
            for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
                for tag in constants.tags:
                    halo = assembily_history[sim][tag]["halo"]
                    subhalo = assembily_history[sim][tag]["subhalo"]
                    if halo is None:
                        continue
                    try:
                        j = sphere_specific_angular_momentum(halo, subhalo, tag, sim, [particle_type], r_200_fraction)
                    except LookupError:
                        continue
                    j_values_tags[particle_type][r_200_fraction][sim].append(tag)
                    j_values[particle_type][r_200_fraction][sim].append(j)

    with open("j_values.txt", "w") as f:
        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
            f.write("--|| STYPE ||-- {}\n".format(Simulations.to_string(sim)))
            for r_200_fraction in (1.0, 0.1):
                f.write(f"--|| RATIO ||-- {r_200_fraction}\n")
                for particle_type in (ParticleType.gas, ParticleType.star, ParticleType.dark_matter):
                    f.write("--|| PTYPE ||-- {}\n".format(ParticleType.to_string(particle_type)))
                    for i in range(len(j_values_tags[particle_type][r_200_fraction][sim])):
                        f.write(f"{j_values_tags[particle_type][r_200_fraction][sim][i]}:{j_values[particle_type][r_200_fraction][sim][i]}\n")
                    



    gas_particle_numbers = { 1.0 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} },
                             0.1 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} }
                            }
    gas_particle_mass = { 1.0 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} },
                          0.1 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} }
                         }
                        
    for r_200_fraction in (1.0, 0.1):
        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
            IDs_by_type = load_particle_IDs(sim, r_200_fraction)
            star_IDs_at_z0 = IDs_by_type[ParticleType.star]
            for tag in constants.tags:
                halo = assembily_history[sim][tag]["halo"]
                subhalo = assembily_history[sim][tag]["subhalo"]
                if halo is None:
                    continue
                try:
                    snapshot = ParticleReadConversion_EagleSnapshot(tag, sim, SimulationModels.RECAL, relitive_data_root)
                    gas_IDs = snapshot.particle_read(ParticleType.gas, "ParticleIDs", unit_system = UnitSystem.cgs)
                    gas_masses = snapshot.particle_read(ParticleType.gas, "Mass", unit_system = UnitSystem.cgs)
                    id_filter = np.in1d(gas_IDs, star_IDs_at_z0)
                    num_as_gas_particles = len(gas_IDs[id_filter])
                    gas_mass = sum(gas_masses[id_filter])
                except LookupError:
                    continue
                gas_particle_numbers[r_200_fraction][sim][tag] = num_as_gas_particles
                gas_particle_mass[r_200_fraction][sim][tag] = gas_mass

    with open("star_forming_gas_values.txt", "w") as f:
        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
            f.write("--|| STYPE ||-- {}\n".format(Simulations.to_string(sim)))
            for r_200_fraction in (1.0, 0.1):
                f.write(f"--|| RATIO ||-- {r_200_fraction}\n")
                for tag in list(gas_particle_numbers[r_200_fraction][sim].keys()):
                    f.write(f"{tag}:{gas_particle_numbers[r_200_fraction][sim][tag]}:{gas_particle_mass[r_200_fraction][sim][tag]}\n")
