import numpy as np
from matplotlib import pyplot as plt
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, baryonic_particle_types, baryonic_low_mass_particle_types, UnitSystem
from regional_specific_angular_momentum import set_R_200_fraction
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs, X_Axis_Value, line_colours, line_styles

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

max_r_tags = {         ParticleType.gas : { 1.0 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None },
                                            0.1 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None } },
                      ParticleType.star : { 1.0 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None },
                                            0.1 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None } },
               ParticleType.dark_matter : { 1.0 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None },
                                            0.1 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None } }
              }

max_j_tags = {         ParticleType.gas : { 1.0 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None },
                                            0.1 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None } },
                      ParticleType.star : { 1.0 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None },
                                            0.1 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None } },
               ParticleType.dark_matter : { 1.0 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None },
                                            0.1 : { Simulations.Early : None, Simulations.Organic : None, Simulations.Late : None } }
              }

with open("r_values.txt", "r") as f:
    current_sim = None
    current_ratio = None
    current_particle = None
    for line in f.readlines():
        line = line.rstrip("\n")
        if line[:5] == "--|| ":
            if line[5:10] == "STYPE":
                current_sim = Simulations.from_string(line[16:])
            elif line[5:10] == "RATIO":
                current_ratio = float(line[16:])
            elif line[5:10] == "PTYPE":
                current_particle = ParticleType.from_string(line[16:])
            else: raise ValueError()
        else:
            if current_sim is None or current_ratio is None or current_particle is None:
                raise TypeError()
            tag, value = line.split(":")
            r_values_tags[current_particle][current_ratio][current_sim].append(tag)
            r_values[current_particle][current_ratio][current_sim].append(float(value))

with open("j_values.txt", "r") as f:
    current_sim = None
    current_ratio = None
    current_particle = None
    for line in f.readlines():
        line = line.rstrip("\n")
        if line[:5] == "--|| ":
            if line[5:10] == "STYPE":
                current_sim = Simulations.from_string(line[16:])
            elif line[5:10] == "RATIO":
                current_ratio = float(line[16:])
            elif line[5:10] == "PTYPE":
                current_particle = ParticleType.from_string(line[16:])
            else: raise ValueError()
        else:
            if current_sim is None or current_ratio is None or current_particle is None:
                raise TypeError()
            tag, value = line.split(":")
            j_values_tags[current_particle][current_ratio][current_sim].append(tag)
            j_values[current_particle][current_ratio][current_sim].append(float(value))

for particle_type in list(r_values_tags.keys()):
    for R_200_ratio in list(r_values_tags[particle_type].keys()):
        for sim in list(r_values_tags[particle_type][R_200_ratio].keys()):
            max_r_tags[particle_type][R_200_ratio][sim] = r_values_tags[particle_type][R_200_ratio][sim][r_values[particle_type][R_200_ratio][sim].index(max(r_values[particle_type][R_200_ratio][sim]))]

for particle_type in list(j_values_tags.keys()):
    for R_200_ratio in list(j_values_tags[particle_type].keys()):
        for sim in list(j_values_tags[particle_type][R_200_ratio].keys()):
            max_j_tags[particle_type][R_200_ratio][sim] = j_values_tags[particle_type][R_200_ratio][sim][j_values[particle_type][R_200_ratio][sim].index(max(j_values[particle_type][R_200_ratio][sim]))]

def get_value_at_tag(values_dict, tags_dict, tag, simulation, particle_type, R_200_fraction):
    try:
        return values_dict[particle_type][R_200_fraction][simulation][tags_dict[particle_type][R_200_fraction][simulation].index(tag)]
    except ValueError:
        raise LookupError("No value at tag.")



gas_particle_numbers = { 1.0 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} },
                         0.1 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} }
                        }
gas_particle_mass = { 1.0 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} },
                      0.1 : { Simulations.Early : {}, Simulations.Organic : {}, Simulations.Late : {} }
                     }

with open("star_forming_gas_values.txt", "r") as f:
    current_sim = None
    current_ratio = None
    for line in f.readlines():
        line = line.rstrip("\n")
        if line[:5] == "--|| ":
            if line[5:10] == "STYPE":
                current_sim = Simulations.from_string(line[16:])
            elif line[5:10] == "RATIO":
                current_ratio = float(line[16:])
            else: raise ValueError()
        else:
            if current_sim is None or current_ratio is None:
                raise TypeError()
            tag, number, mass = line.split(":")
            gas_particle_numbers[current_ratio][current_sim][tag] = int(number)
            gas_particle_mass[current_ratio][current_sim][tag] = float(mass)



def set_particle_type(func, particle_type):
    def wrapper(*args, **kwargs):
        return func(particle_type = particle_type, *args, **kwargs)
    return wrapper

def sphere_specific_angular_momentum(halo, subhalo, tag, simulation, particle_type, R_200_fraction = 1.0):
    if tag not in j_values_tags[particle_type][R_200_fraction][simulation]:
        raise LookupError("Value not in cache.")
    return get_value_at_tag(j_values, j_values_tags, tag, simulation, particle_type, R_200_fraction)

def sphere_radius(halo, subhalo, tag, simulation, particle_type, R_200_fraction = 1.0):
    if tag not in r_values_tags[particle_type][R_200_fraction][simulation]:
        raise LookupError("Value not in cache.")
    return get_value_at_tag(r_values, r_values_tags, tag, simulation, particle_type, R_200_fraction)

def j_ratio(halo, subhalo, tag, simulation, particle_type, R_200_fraction = 1.0):
    j_max_tag = max_j_tags[particle_type][R_200_fraction][simulation]
    j_max = get_value_at_tag(j_values, j_values_tags, j_max_tag, simulation, particle_type, R_200_fraction)
    j_now = get_value_at_tag(j_values, j_values_tags, tag, simulation, particle_type, R_200_fraction)

    return j_now / j_max

def star_formaing_gas_fraction(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    earliest_tag = None
    for earliest_tag in constants.tags:
        if earliest_tag in list(gas_particle_mass[R_200_fraction][simulation].keys()):
            break
    return gas_particle_mass[R_200_fraction][simulation][tag] / gas_particle_mass[R_200_fraction][simulation][earliest_tag]



if __name__ == "__main__":
    # j Dark Matter with Markers
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.dark_matter][1.0][simulation]
            j_at_r_max = get_value_at_tag(j_values, j_values_tags, r_max_tag, simulation, ParticleType.dark_matter, 1.0)
            plt.plot([constants.expansion_factors[r_max_tag]], [j_at_r_max], marker = "s", color = line_colours[3], label = "Max R" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.dark_matter][1.0][simulation]
            j_max = get_value_at_tag(j_values, j_values_tags, j_max_tag, simulation, ParticleType.dark_matter, 1.0)
            plt.plot([constants.expansion_factors[j_max_tag]], [j_max], marker = "o", color = line_colours[3], label = "Max J" if i == 0 else None)
    produce_simulations_graph(set_particle_type(sphere_specific_angular_momentum, ParticleType.dark_matter),
                              "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = plot_maxes,
                              filename = "j_DM_R200_with_markers.png")
    
    # j Dark Matter with Markers (inner halo)
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.dark_matter][0.1][simulation]
            j_at_r_max = get_value_at_tag(j_values, j_values_tags, r_max_tag, simulation, ParticleType.dark_matter, 0.1)
            plt.plot([constants.expansion_factors[r_max_tag]], [j_at_r_max], marker = "s", color = line_colours[3], label = "Max R" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.dark_matter][0.1][simulation]
            j_max = get_value_at_tag(j_values, j_values_tags, j_max_tag, simulation, ParticleType.dark_matter, 0.1)
            plt.plot([constants.expansion_factors[j_max_tag]], [j_max], marker = "o", color = line_colours[3], label = "Max J" if i == 0 else None)
    produce_simulations_graph(set_R_200_fraction(set_particle_type(sphere_specific_angular_momentum, ParticleType.dark_matter), 0.1),
                              "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = plot_maxes,
                              filename = "j_DM_0.1_R200_with_markers.png")

    # j_star ratio vs. star as gas fraction
    produce_single_simulation_graphs([set_particle_type(j_ratio, ParticleType.star), star_formaing_gas_fraction],
                                     ["Fraction of Maximum |$\\vec{j_{star}}$|", "Fraction of Unconverted Gas Mass"], y_axis_label = "Fraction", graph_title_partial = "Fraction of Maximum Achieved Stellar Specific Angular Momentum (particles within $R_{200}$ at z=0)",
                                     x_axis = X_Axis_Value.time, log_x = False, log_y = False, invert_x = False, use_rolling_average = False,
                                     xlim_overide = (list(constants.times.values())[0], list(constants.times.values())[-1]), ylim_overide = (0, 1),
                                     filename = "j_star_R200_fraction_of_max_with_unconverted_gas.png")

    # j_star ratio vs. star as gas fraction (inner halo)
    produce_single_simulation_graphs([set_R_200_fraction(set_particle_type(j_ratio, ParticleType.star), 0.1), set_R_200_fraction(star_formaing_gas_fraction, 0.1)],
                                     ["Fraction of Maximum |$\\vec{j_{star}}$|", "Fraction of Unconverted Gas Mass"], y_axis_label = "Fraction", graph_title_partial = "Fraction of Maximum Achieved Stellar Specific Angular Momentum (particles within 0.1 $R_{200}$ at z=0)",
                                     x_axis = X_Axis_Value.time, log_x = False, log_y = False, invert_x = False, use_rolling_average = False,
                                     xlim_overide = (list(constants.times.values())[0], list(constants.times.values())[-1]), ylim_overide = (0, 1),
                                     filename = "j_star_0.1_R200_fraction_of_max_with_unconverted_gas.png")
