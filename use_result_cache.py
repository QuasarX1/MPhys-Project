import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, baryonic_particle_types, baryonic_low_mass_particle_types, UnitSystem
from regional_specific_angular_momentum import set_R_200_fraction
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs, produce_overlayed_single_simulation_graphs, X_Axis_Value, line_colours, line_styles, plotting_font_size, plotting_line_size

marker_size = 10

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

def j_star(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    return get_value_at_tag(j_values, j_values_tags, tag, simulation, ParticleType.star, R_200_fraction)

def j_DM(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    return get_value_at_tag(j_values, j_values_tags, tag, simulation, ParticleType.dark_matter, R_200_fraction)

def r_star(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    return get_value_at_tag(r_values, r_values_tags, tag, simulation, ParticleType.star, R_200_fraction)

def r_DM(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    return get_value_at_tag(r_values, r_values_tags, tag, simulation, ParticleType.dark_matter, R_200_fraction)

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
    # Specific Angular Momentum Ratio
    print("Stars:")
    print("R_200:\n      Early: {:.4f}\n    Organic: {:.4f}\n       Late: {:.4f}".format(
        j_ratio(None, None, constants.tags[-1], Simulations.Early, ParticleType.star, 1.0),
        j_ratio(None, None, constants.tags[-1], Simulations.Organic, ParticleType.star, 1.0),
        j_ratio(None, None, constants.tags[-1], Simulations.Late, ParticleType.star, 1.0)))
    print("0.1 R_200:\n      Early: {:.4f}\n    Organic: {:.4f}\n       Late: {:.4f}".format(
        j_ratio(None, None, constants.tags[-1], Simulations.Early, ParticleType.star, 0.1),
        j_ratio(None, None, constants.tags[-1], Simulations.Organic, ParticleType.star, 0.1),
        j_ratio(None, None, constants.tags[-1], Simulations.Late, ParticleType.star, 0.1)))
    print("DM:")
    print("R_200:\n      Early: {:.4f}\n    Organic: {:.4f}\n       Late: {:.4f}".format(
        j_ratio(None, None, constants.tags[-1], Simulations.Early, ParticleType.dark_matter, 1.0),
        j_ratio(None, None, constants.tags[-1], Simulations.Organic, ParticleType.dark_matter, 1.0),
        j_ratio(None, None, constants.tags[-1], Simulations.Late, ParticleType.dark_matter, 1.0)))
    print("0.1 R_200:\n      Early: {:.4f}\n    Organic: {:.4f}\n       Late: {:.4f}".format(
        j_ratio(None, None, constants.tags[-1], Simulations.Early, ParticleType.dark_matter, 0.1),
        j_ratio(None, None, constants.tags[-1], Simulations.Organic, ParticleType.dark_matter, 0.1),
        j_ratio(None, None, constants.tags[-1], Simulations.Late, ParticleType.dark_matter, 0.1)))
    exit()

    # j_star ratio vs. star as gas fraction
    def max_point_marker(sim, axis, x_values, y_values, *args, **kwargs):
        x_values = [list(x_values[i]) for i in range(len(x_values))]
        y_values = [list(y_values[i]) for i in range(len(y_values))]

        j_ratio_max_x = x_values[0][y_values[0].index(max(y_values[0]))]
        gas_fraction_y = y_values[1][x_values[1].index(x_values[0][y_values[0].index(max(y_values[0]))])]
        
        axis.plot([j_ratio_max_x, j_ratio_max_x], [0.5, gas_fraction_y], color = line_colours[5])
        axis.plot([j_ratio_max_x], [gas_fraction_y], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "t at max |$\\vec{j_{star}}$|")
        axis.plot([j_ratio_max_x], [0.5], linestyle = "", marker = "s", markersize = marker_size, color = line_colours[4], label = "$\\frac{1}{2}$ $M_{gas}$")
    produce_single_simulation_graphs([set_particle_type(j_ratio, ParticleType.star), star_formaing_gas_fraction],
                                     ["Fraction of Maximum |$\\vec{j_{star}}$|", "Fraction of Unconverted Gas Mass"], y_axis_label = "Fraction", graph_title_partial = "Fraction of Maximum Achieved Stellar\nSpecific Angular Momentum (particles within $R_{200}$ at z=0)",
                                     x_axis = X_Axis_Value.time, log_x = False, log_y = False, invert_x = False, use_rolling_average = False,
                                     extra_plotting_funcs = [lambda *args, **kwargs: max_point_marker(Simulations.Early, *args, **kwargs), lambda *args, **kwargs: max_point_marker(Simulations.Organic, *args, **kwargs), lambda *args, **kwargs: max_point_marker(Simulations.Late, *args, **kwargs)],
                                     xlim_overide = (list(constants.times.values())[0], list(constants.times.values())[-1]), ylim_overide = (0, 1),
                                     filename = "j_star_R200_fraction_of_max_with_unconverted_gas.png", vertical_stack = True)
    #produce_overlayed_single_simulation_graphs([set_particle_type(j_ratio, ParticleType.star), star_formaing_gas_fraction],
    #                                           ["Fraction of Maximum |$\\vec{j_{star}}$|", "Fraction of Unconverted Gas Mass"], y_axis_label = "Fraction", graph_title_partial = "Fraction of Maximum Achieved Stellar Specific Angular Momentum (particles within $R_{200}$ at z=0)",
    #                                           x_axis = X_Axis_Value.time, log_x = False, log_y = False, invert_x = False, use_rolling_average = False,
    #                                           xlim_overide = (list(constants.times.values())[0], list(constants.times.values())[-1]), ylim_overide = (0, 1),
    #                                           filename = "j_star_R200_fraction_of_max_with_unconverted_gas_overlapped.png")

    # j_star ratio vs. star as gas fraction (inner halo)
    def max_point_marker(sim, axis, x_values, y_values, *args, **kwargs):
        x_values = [list(x_values[i]) for i in range(len(x_values))]
        y_values = [list(y_values[i]) for i in range(len(y_values))]

        j_ratio_max_x = x_values[0][y_values[0].index(max(y_values[0]))]
        gas_fraction_y = y_values[1][x_values[1].index(x_values[0][y_values[0].index(max(y_values[0]))])]
        
        axis.plot([j_ratio_max_x, j_ratio_max_x], [0.5, gas_fraction_y], color = line_colours[5])
        axis.plot([j_ratio_max_x], [gas_fraction_y], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "t at max |$\\vec{j_{star}}$|")
        axis.plot([j_ratio_max_x], [0.5], linestyle = "", marker = "s", markersize = marker_size, color = line_colours[4], label = "$\\frac{1}{2}$ $M_{gas}$")
    produce_single_simulation_graphs([set_R_200_fraction(set_particle_type(j_ratio, ParticleType.star), 0.1), set_R_200_fraction(star_formaing_gas_fraction, 0.1)],
                                     ["Fraction of Maximum |$\\vec{j_{star}}$|", "Fraction of Unconverted Gas Mass"], y_axis_label = "Fraction", graph_title_partial = "Fraction of Maximum Achieved Stellar\nSpecific Angular Momentum (particles within 0.1 $R_{200}$ at z=0)",
                                     x_axis = X_Axis_Value.time, log_x = False, log_y = False, invert_x = False, use_rolling_average = False,
                                     extra_plotting_funcs = [lambda *args, **kwargs: max_point_marker(Simulations.Early, *args, **kwargs), lambda *args, **kwargs: max_point_marker(Simulations.Organic, *args, **kwargs), lambda *args, **kwargs: max_point_marker(Simulations.Late, *args, **kwargs)],
                                     xlim_overide = (list(constants.times.values())[0], list(constants.times.values())[-1]), ylim_overide = (0, 1),
                                     filename = "j_star_0.1_R200_fraction_of_max_with_unconverted_gas.png", vertical_stack = True)
    #produce_overlayed_single_simulation_graphs([set_R_200_fraction(set_particle_type(j_ratio, ParticleType.star), 0.1), set_R_200_fraction(star_formaing_gas_fraction, 0.1)],
    #                                           ["Fraction of Maximum |$\\vec{j_{star}}$|", "Fraction of Unconverted Gas Mass"], y_axis_label = "Fraction", graph_title_partial = "Fraction of Maximum Achieved Stellar Specific Angular Momentum (particles within 0.1 $R_{200}$ at z=0)",
    #                                           x_axis = X_Axis_Value.time, log_x = False, log_y = False, invert_x = False, use_rolling_average = False,
    #                                           xlim_overide = (list(constants.times.values())[0], list(constants.times.values())[-1]), ylim_overide = (0, 1),
    #                                           filename = "j_star_0.1_R200_fraction_of_max_with_unconverted_gas_overlapped.png")



    # Specific Angular Momentum
    # j_star
    produce_simulations_graph(set_R_200_fraction(j_star, 1.0),
                              y_axis_label = "|$\\vec{j_{star}}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum of Star particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size),
                              filename = "j_star_R200.png")
    # With Markers
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            j_DM_max_tag = max_j_tags[ParticleType.dark_matter][1.0][simulation]
            j_at_max_DM = get_value_at_tag(j_values, j_values_tags, j_DM_max_tag, simulation, ParticleType.star, 1.0)
            #dy = 3 * 10**(np.log10(j_at_max_DM) - 1)
            #plt.arrow(constants.expansion_factors[j_DM_max_tag], j_at_max_DM - dy, 0, dy, width = 0.005, label = "Max |$\\vec{j_{DM}}$|" if i == 0 else None)
            plt.plot([constants.expansion_factors[j_DM_max_tag]], [j_at_max_DM], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j_{DM}}$|" if i == 0 else None)
            #r_max_tag = max_r_tags[ParticleType.star][1.0][simulation]
            #j_at_r_max = get_value_at_tag(j_values, j_values_tags, r_max_tag, simulation, ParticleType.star, 1.0)
            #plt.plot([constants.expansion_factors[r_max_tag]], [j_at_r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[-1], label = "a at max median($r_{star}$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.star][1.0][simulation]
            j_max = get_value_at_tag(j_values, j_values_tags, j_max_tag, simulation, ParticleType.star, 1.0)
            plt.plot([constants.expansion_factors[j_max_tag]], [j_max], linestyle = "", marker = "*", markersize = marker_size, color = line_colours[-1], label = "a at max |$\\vec{j_{star}}$|" if i == 0 else None)
        plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size)
    produce_simulations_graph(set_R_200_fraction(j_star, 1.0),
                              y_axis_label = "|$\\vec{j_{star}}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum of Star particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = plot_maxes,
                              filename = "j_star_R200_with_markers.png")

    # j_star inner
    produce_simulations_graph(set_R_200_fraction(j_star, 0.1),
                              y_axis_label = "|$\\vec{j_{star}}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum of Star particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size),
                              filename = "j_star_0.1_R200.png")
    # With Markers (inner halo)
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            j_DM_max_tag = max_j_tags[ParticleType.dark_matter][0.1][simulation]
            j_at_max_DM = get_value_at_tag(j_values, j_values_tags, j_DM_max_tag, simulation, ParticleType.star, 0.1)
            #dy = 3 * 10**(np.log10(j_at_max_DM) - 1)
            #plt.arrow(constants.expansion_factors[j_DM_max_tag], j_at_max_DM - dy, 0, dy)
            plt.plot([constants.expansion_factors[j_DM_max_tag]], [j_at_max_DM], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j_{DM}}$|" if i == 0 else None)
            #r_max_tag = max_r_tags[ParticleType.star][0.1][simulation]
            #j_at_r_max = get_value_at_tag(j_values, j_values_tags, r_max_tag, simulation, ParticleType.star, 0.1)
            #plt.plot([constants.expansion_factors[r_max_tag]], [j_at_r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[-1], label = "a at max median($r_{star}$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.star][0.1][simulation]
            j_max = get_value_at_tag(j_values, j_values_tags, j_max_tag, simulation, ParticleType.star, 0.1)
            plt.plot([constants.expansion_factors[j_max_tag]], [j_max], linestyle = "", marker = "*", markersize = marker_size, color = line_colours[-1], label = "a at max |$\\vec{j_{star}}$|" if i == 0 else None)
        plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size)
    produce_simulations_graph(set_R_200_fraction(j_star, 0.1),
                              y_axis_label = "|$\\vec{j_{star}}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum of Star particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = plot_maxes,
                              filename = "j_star_0.1_R200_with_markers.png")

    # Dark Matter
    produce_simulations_graph(set_R_200_fraction(j_DM, 1.0),
                              y_axis_label = "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum of Dark Matter particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size),
                              filename = "j_DM_R200.png")
    # With Markers
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.dark_matter][1.0][simulation]
            j_at_r_max = get_value_at_tag(j_values, j_values_tags, r_max_tag, simulation, ParticleType.dark_matter, 1.0)
            plt.plot([constants.expansion_factors[r_max_tag]], [j_at_r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[3], label = "a at max median($r$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.dark_matter][1.0][simulation]
            j_max = get_value_at_tag(j_values, j_values_tags, j_max_tag, simulation, ParticleType.dark_matter, 1.0)
            plt.plot([constants.expansion_factors[j_max_tag]], [j_max], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j}$|" if i == 0 else None)
        plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size)
    produce_simulations_graph(set_particle_type(sphere_specific_angular_momentum, ParticleType.dark_matter),
                              "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = plot_maxes,
                              filename = "j_DM_R200_with_markers.png")

    # Dark Matter inner
    produce_simulations_graph(set_R_200_fraction(j_DM, 0.1),
                              y_axis_label = "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum of Dark Matter particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size),
                              filename = "j_DM_0.1_R200.png")
    # With Markers (inner halo)
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.dark_matter][0.1][simulation]
            j_at_r_max = get_value_at_tag(j_values, j_values_tags, r_max_tag, simulation, ParticleType.dark_matter, 0.1)
            plt.plot([constants.expansion_factors[r_max_tag]], [j_at_r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[3], label = "a at max median($r$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.dark_matter][0.1][simulation]
            j_max = get_value_at_tag(j_values, j_values_tags, j_max_tag, simulation, ParticleType.dark_matter, 0.1)
            plt.plot([constants.expansion_factors[j_max_tag]], [j_max], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j}$|" if i == 0 else None)
        plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[3], linestyle = line_styles[-1], linewidth = plotting_line_size)
    produce_simulations_graph(set_R_200_fraction(set_particle_type(sphere_specific_angular_momentum, ParticleType.dark_matter), 0.1),
                              "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                              extra_plotting_func = plot_maxes,
                              filename = "j_DM_0.1_R200_with_markers.png")



    # Median Radius
    # r_star
    produce_simulations_graph(set_R_200_fraction(r_star, 1.0),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Star particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (10, 500),
                              filename = "r_star_R200.png")
    # With Markers
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.star][1.0][simulation]
            r_max = get_value_at_tag(r_values, r_values_tags, r_max_tag, simulation, ParticleType.star, 1.0)
            plt.plot([constants.expansion_factors[r_max_tag]], [r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[3], label = "a at max median($r$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.star][1.0][simulation]
            r_at_j_max = get_value_at_tag(r_values, r_values_tags, j_max_tag, simulation, ParticleType.star, 1.0)
            plt.plot([constants.expansion_factors[j_max_tag]], [r_at_j_max], linestyle = "", marker = "*", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j}$|" if i == 0 else None)
    produce_simulations_graph(set_R_200_fraction(r_star, 1.0),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Star particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (10, 500),
                              extra_plotting_func = plot_maxes,
                              filename = "r_star_R200_with_markers.png")

    # r_star inner
    produce_simulations_graph(set_R_200_fraction(r_star, 0.1),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Star particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (10, 500),
                              filename = "r_star_0.1_R200.png")
    # With Markers (inner halo)
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.star][0.1][simulation]
            r_max = get_value_at_tag(r_values, r_values_tags, r_max_tag, simulation, ParticleType.star, 0.1)
            plt.plot([constants.expansion_factors[r_max_tag]], [r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[3], label = "a at max median($r$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.star][0.1][simulation]
            r_at_j_max = get_value_at_tag(r_values, r_values_tags, j_max_tag, simulation, ParticleType.star, 0.1)
            plt.plot([constants.expansion_factors[j_max_tag]], [r_at_j_max], linestyle = "", marker = "*", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j}$|" if i == 0 else None)
    produce_simulations_graph(set_R_200_fraction(r_star, 0.1),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Star particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (10, 500),
                              extra_plotting_func = plot_maxes,
                              filename = "r_star_0.1_R200_with_markers.png")

    # r_DM
    produce_simulations_graph(set_R_200_fraction(r_DM, 1.0),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Dark Matter particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (100, 500),
                              extra_plotting_func = lambda x_values, *args, **kwargs: plt.yticks([100, 200, 300, 400, 500], ["100", "200", "300", "400", "500"], fontsize = plotting_font_size),
                              filename = "r_DM_R200.png")
    # With Markers
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.dark_matter][1.0][simulation]
            r_max = get_value_at_tag(r_values, r_values_tags, r_max_tag, simulation, ParticleType.dark_matter, 1.0)
            plt.plot([constants.expansion_factors[r_max_tag]], [r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[3], label = "a at max median($r$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.dark_matter][1.0][simulation]
            r_at_j_max = get_value_at_tag(r_values, r_values_tags, j_max_tag, simulation, ParticleType.dark_matter, 1.0)
            plt.plot([constants.expansion_factors[j_max_tag]], [r_at_j_max], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j}$|" if i == 0 else None)
        plt.yticks([100, 200, 300, 400, 500], ["$10^2$", "", "3$\\times10^2$", "", "5$\\times10^2$"], fontsize = plotting_font_size)
        #plt.yticks([100, 300, 500], ["100", "3 $\\times$ $10^2$", "5 $\\times$ $10^2$"], fontsize = plotting_font_size)
        #plt.gca().get_yaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    produce_simulations_graph(set_R_200_fraction(r_DM, 1.0),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Dark Matter particles within $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              #xlim_overide = (0.1, 1), ylim_overide = (10, 500),
                              xlim_overide = (0.1, 1), ylim_overide = (100, 500),
                              extra_plotting_func = plot_maxes,
                              filename = "r_DM_R200_with_markers.png")
        
    # r_DM inner
    produce_simulations_graph(set_R_200_fraction(r_DM, 0.1),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Dark Matter particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (10, 500),
                              filename = "r_DM_0.1_R200.png")
    # With Markers (inner halo)
    def plot_maxes(x_values, *args, **kwargs):
        for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
            r_max_tag = max_r_tags[ParticleType.dark_matter][0.1][simulation]
            r_max = get_value_at_tag(r_values, r_values_tags, r_max_tag, simulation, ParticleType.dark_matter, 0.1)
            plt.plot([constants.expansion_factors[r_max_tag]], [r_max], linestyle = "", marker = "^", markersize = marker_size, color = line_colours[3], label = "a at max median($r$)" if i == 0 else None)
            j_max_tag = max_j_tags[ParticleType.dark_matter][0.1][simulation]
            r_at_j_max = get_value_at_tag(r_values, r_values_tags, j_max_tag, simulation, ParticleType.dark_matter, 0.1)
            plt.plot([constants.expansion_factors[j_max_tag]], [r_at_j_max], linestyle = "", marker = "o", markersize = marker_size, color = line_colours[3], label = "a at max |$\\vec{j}$|" if i == 0 else None)
    produce_simulations_graph(set_R_200_fraction(r_DM, 0.1),
                              y_axis_label = "median($r$) ($kPc$)", graph_title_partial = "Median Radius of Dark Matter particles within 0.1 $R_{200}$ at z=0",
                              x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                              xlim_overide = (0.1, 1), ylim_overide = (10, 500),
                              extra_plotting_func = plot_maxes,
                              filename = "r_DM_0.1_R200_with_markers.png")
        