import numpy as np
from matplotlib import pyplot as plt
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, baryonic_particle_types, baryonic_low_mass_particle_types, UnitSystem
from assembily_history import physical_centeral_mass_positions, assembily_history
from z0_particle_IDs import load_particle_IDs
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs, X_Axis_Value, line_colours, line_styles
from Physics import specific_angular_momentum, centre_of_mass
from angular_momentum_units import specific_angular_momentum_units

relitive_data_root = ".\\gm_for_mphys"

def set_particle_types(func, *particle_types):
    def wrapper(*args, **kwargs):
        return func(particle_types = particle_types, *args, **kwargs)
    return wrapper

def set_R_200_fraction(func, fraction):
    def wrapper(*args, **kwargs):
        return func(R_200_fraction = fraction, *args, **kwargs)
    return wrapper

def sphere_specific_angular_momentum(halo, subhalo, tag, simulation, particle_types = all_particle_types, R_200_fraction = 1.0):
    print("Doing {} {}".format(tag, Simulations.to_string(simulation)))
    #galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    #if galaxy_centre is None:
    #    raise LookupError("No data avalible.")


    particle_IDs_by_type = load_particle_IDs(simulation, R_200_fraction)
    # Use below for using the same particles as organic case in all cases
    #particle_IDs_by_type = load_particle_IDs(Simulations.Organic, R_200_fraction)


    relivant_particle_IDs = np.empty((sum([len(particle_IDs_by_type[particle_type]) for particle_type in particle_types]), ))
    next_index = 0
    for particle_type in particle_types:
        relivant_particle_IDs[next_index : next_index + len(particle_IDs_by_type[particle_type])] = particle_IDs_by_type[particle_type]
        next_index += len(particle_IDs_by_type[particle_type])

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    box_all_positions = np.empty(shape = (0, 3))
    box_all_velocities = np.empty(shape = (0, 3))
    box_all_masses = np.empty(shape = (0, ))
    for particle_type in particle_IDs_by_type:
        box_particle_IDs, box_particle_locations = snapshot.particle_read(particle_type, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)
        id_filter = np.in1d(box_particle_IDs, relivant_particle_IDs)
        if sum(id_filter) == 0:
            continue# No particles with matching IDs present in this particle type
        box_particle_locations = box_particle_locations[id_filter]
        box_particle_velocities = snapshot.particle_read(particle_type, "Velocity", unit_system = UnitSystem.cgs)[id_filter]
        if particle_type != ParticleType.dark_matter:
            box_particle_masses = snapshot.particle_read(particle_type, "Mass", unit_system = UnitSystem.cgs)[id_filter]
        else:
            box_particle_masses = np.full(box_particle_locations.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        box_all_positions = np.append(box_all_positions, box_particle_locations, axis = 0)
        box_all_velocities = np.append(box_all_velocities, box_particle_velocities, axis = 0)
        box_all_masses = np.append(box_all_masses, box_particle_masses)
        
    selection_centre = centre_of_mass(box_all_positions, box_all_masses)
    all_positions = snapshot.convert_distance_values(
                        snapshot.centre_particles(
                            snapshot.convert_distance_values(box_all_positions, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET),
                            snapshot.convert_distance_values(selection_centre, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET)
                        ),
                        UnitSystem.h_less_comoving_GADGET, UnitSystem.cgs
                    )

    return specific_angular_momentum_units(np.sqrt((specific_angular_momentum(all_positions, box_all_velocities, box_all_masses)**2).sum()))

def sphere_radius(halo, subhalo, tag, simulation, particle_types = all_particle_types, R_200_fraction = 1.0):
    print("Doing {} {}".format(tag, Simulations.to_string(simulation)))
    #galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    #if galaxy_centre is None:
    #    raise LookupError("No data avalible.")


    particle_IDs_by_type = load_particle_IDs(simulation, R_200_fraction)
    # Use below for using the same particles as organic case in all cases
    #particle_IDs_by_type = load_particle_IDs(Simulations.Organic, R_200_fraction)


    relivant_particle_IDs = np.empty((sum([len(particle_IDs_by_type[particle_type]) for particle_type in particle_types]), ))
    next_index = 0
    for particle_type in particle_types:
        relivant_particle_IDs[next_index : next_index + len(particle_IDs_by_type[particle_type])] = particle_IDs_by_type[particle_type]
        next_index += len(particle_IDs_by_type[particle_type])

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    box_all_positions = np.empty(shape = (0, 3))
    box_all_masses = np.empty(shape = (0, ))
    for particle_type in particle_IDs_by_type:
        box_particle_IDs, box_particle_locations = snapshot.particle_read(particle_type, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)
        id_filter = np.in1d(box_particle_IDs, relivant_particle_IDs)
        if sum(id_filter) == 0:
            continue# No particles with matching IDs present in this particle type
        box_particle_locations = box_particle_locations[id_filter]
        if particle_type != ParticleType.dark_matter:
            box_particle_masses = snapshot.particle_read(particle_type, "Mass", unit_system = UnitSystem.cgs)[id_filter]
        else:
            box_particle_masses = np.full(box_particle_locations.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        box_all_positions = np.append(box_all_positions, box_particle_locations, axis = 0)
        box_all_masses = np.append(box_all_masses, box_particle_masses)
        
    selection_centre = centre_of_mass(box_all_positions, box_all_masses)
    all_positions = snapshot.convert_distance_values(
                        snapshot.centre_particles(
                            snapshot.convert_distance_values(box_all_positions, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET),
                            snapshot.convert_distance_values(selection_centre, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET)
                        ),
                        UnitSystem.h_less_comoving_GADGET, UnitSystem.cgs
                    )

    radi = np.sqrt(np.sum(all_positions**2, axis = 1))
    
    return np.median(radi) * 10**(3) / constants.SimulationConstants.get_constants()["CM_PER_MPC"]

#particle_IDs_by_type = load_particle_IDs(Simulations.Early, 1.0)
#print("Early:")
#print("N Gas Particles: " + str(len(particle_IDs_by_type[ParticleType.gas])))
#print("N Star Particles: " + str(len(particle_IDs_by_type[ParticleType.star])))
#particle_IDs_by_type = load_particle_IDs(Simulations.Organic, 1.0)
#print("Organic:")
#print("N Gas Particles: " + str(len(particle_IDs_by_type[ParticleType.gas])))
#print("N Star Particles: " + str(len(particle_IDs_by_type[ParticleType.star])))
#particle_IDs_by_type = load_particle_IDs(Simulations.Late, 1.0)
#print("Late:")
#print("N Gas Particles: " + str(len(particle_IDs_by_type[ParticleType.gas])))
#print("N Star Particles: " + str(len(particle_IDs_by_type[ParticleType.star])))
#exit()

if __name__ == "__main__":
    # 1.0 R_200 ----------------------------------------------------------------

    ## Dark Matter
    #print("Dark Matter - j")
    #produce_simulations_graph(set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter),
    #                          "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]),
    #                          filename = "j_DM_R200.png")
    #print("Dark Matter - r")
    #produce_simulations_graph(set_particle_types(sphere_radius, ParticleType.dark_matter),
    #                          "median($r$) ($kPc$)", "Median Radius of Dark Matter particles within $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (10, 500),
    #                          filename = "r_all_DM_R200.png")

    ## Baryons (no Black Holes)
    #print("Baryons - j")
    #produce_simulations_graph(set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types),
    #                          "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Baryon Particles within $R_{200}$ at z=0 (excluding Black Holes)",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]),
    #                          filename = "j_baryon_R200.png")
    #print("Gas - r")
    #produce_simulations_graph(set_particle_types(sphere_radius, ParticleType.gas),
    #                          "median($r$) ($kPc$)", "Median Radius of Gas particles within $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (10, 500),
    #                          filename = "r_gas_R200.png")
    #print("Stars - r")
    #produce_simulations_graph(set_particle_types(sphere_radius, ParticleType.star),
    #                          "median($r$) ($kPc$)", "Median Radius of Star particles within $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (10, 500),
    #                          filename = "r_star_R200.png")

    ## Dark Matter vs. Baryons
    #print("Dark Matter vs. Baryons - j")
    #produce_single_simulation_graphs([set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types)],
    #                                 ["Dark Matter", "Baryonic Matter"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within $R_{200}$ at z=0)",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_DM_to_baryon_R200_comparison.png")

    ## Dark Matter vs. Baryons vs. Stars
    #print("Dark Matter vs. Baryons vs. Stars - j")
    #produce_single_simulation_graphs([set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types), set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter.star)],
    #                                 ["Dark Matter", "Baryonic Matter", "Stars"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within $R_{200}$ at z=0)",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_DM_to_baryon_to_star_R200_comparison.png")

    ## Dark Matter vs. Gas vs. Stars
    #print("Dark Matter vs. Baryons vs. Stars - j")
    #produce_single_simulation_graphs([set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), set_particle_types(sphere_specific_angular_momentum, ParticleType.gas), set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter.star)],
    #                                 ["Dark Matter", "Gas", "Stars"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within $R_{200}$ at z=0)",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_DM_to_gas_to_star_R200_comparison.png")

    ## Baryons / Dark Matter
    #print("Dark Matter / Baryons - j")
    #produce_simulations_graph(lambda *args, **kwargs: sphere_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, *args, **kwargs) / sphere_specific_angular_momentum(particle_types = [ParticleType.dark_matter], *args, **kwargs),
    #                          "|$\\vec{j_b}$ / $\\vec{j_DM}$| (Arbitrary Units)", "Ratio of Baryonic to Dark Matter Specific Angular Momentum for particles within $R_{200}$ at z=0 (excluding Black Holes)",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = False, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1),
    #                          filename = "j_baryon_to_DM_fraction_R200.png")



    # 0.1 R_200 ----------------------------------------------------------------

    ## Dark Matter
    #print("Dark Matter - j - 0.1 R_200")
    #produce_simulations_graph(set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), 0.1),
    #                          "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within 0.1 $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]),
    #                          filename = "j_DM_0.1_R200.png")
    #print("Dark Matter - r - 0.1 R_200")
    #produce_simulations_graph(set_R_200_fraction(set_particle_types(sphere_radius, ParticleType.dark_matter), 0.1),
    #                          "median($r$) ($kPc$)", "Median Radius of Dark Matter particles within 0.1 $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (10, 500),
    #                          filename = "r_all_DM_0.1_R200.png")

    ## Baryons (no Black Holes)
    #print("Baryons - j - 0.1 R_200")
    #produce_simulations_graph(set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types), 0.1),
    #                          "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Baryon Particles within 0.1 $R_{200}$ at z=0 (excluding Black Holes)",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $a^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]),
    #                          filename = "j_baryon_0.1_R200.png")
    #print("Gas - r - 0.1 R_200")
    #produce_simulations_graph(set_R_200_fraction(set_particle_types(sphere_radius, ParticleType.gas), 0.1),
    #                          "median($r$) ($kPc$)", "Median Radius of Gas particles within 0.1 $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (10, 500),
    #                          filename = "r_gas_0.1_R200.png")
    #print("Stars - r - 0.1 R_200")
    #produce_simulations_graph(set_R_200_fraction(set_particle_types(sphere_radius, ParticleType.star), 0.1),
    #                          "median($r$) ($kPc$)", "Median Radius of Star particles within 0.1 $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1), ylim_overide = (10, 500),
    #                          filename = "r_star_0.1_R200.png")

    ## Dark Matter vs. Baryons
    #print("Dark Matter vs. Baryons - j - 0.1 R_200")
    #produce_single_simulation_graphs([set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), 0.1), set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types), 0.1)],
    #                                 ["Dark Matter", "Baryonic Matter"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within 0.1 $R_{200}$ at z=0)",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_DM_to_baryon_0.1_R200_comparison.png")

    ## Dark Matter vs. Baryons vs. Stars
    #print("Dark Matter vs. Baryons vs. Stars - j - 0.1 R_200")
    #produce_single_simulation_graphs([set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), 0.1), set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types), 0.1), set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.star), 0.1)],
    #                                 ["Dark Matter", "Baryonic Matter", "Stars"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within 0.1 $R_{200}$ at z=0)",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_DM_to_baryon_to_star_0.1_R200_comparison.png")

    ## Dark Matter vs. Gas vs. Stars
    #print("Dark Matter vs. Gas vs. Stars - j - 0.1 R_200")
    #produce_single_simulation_graphs([set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), 0.1), set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.gas), 0.1), set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.star), 0.1)],
    #                                 ["Dark Matter", "Gas", "Stars"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within 0.1 $R_{200}$ at z=0)",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_DM_to_gas_to_star_0.1_R200_comparison.png")

    ## Baryons / Dark Matter
    #print("Dark Matter / Baryons - j - 0.1 R_200")
    #produce_simulations_graph(lambda *args, **kwargs: sphere_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, R_200_fraction = 0.1, *args, **kwargs) / sphere_specific_angular_momentum(particle_types = [ParticleType.dark_matter], R_200_fraction = 0.1, *args, **kwargs),
    #                          "|$\\vec{j_b}$ / $\\vec{j_DM}$| (Arbitrary Units)", "Ratio of Baryonic to Dark Matter Specific Angular Momentum for particles within 0.1 $R_{200}$ at z=0 (excluding Black Holes)",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = False, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1),
    #                          filename = "j_baryon_to_DM_fraction_0.1_R200.png")



    # 0.1 R_200 / 1.0 R_200 ----------------------------------------------------

    ## Dark Matter Comparison
    #print("Dark Matter - j - 0.1 R_200 and R_200")
    #produce_single_simulation_graphs([set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), 1.0), set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, ParticleType.dark_matter), 0.1)],
    #                                 ["R = $R_{200}$", "R = 0.1 $R_{200}$"], y_axis_label = "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Dark Matter Specific Angular Momentum by particles within radius R at z=0",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_DM_inner_halo_and_halo_comparison.png")

    ## Baryons Comparison
    #print("Baryons - j - 0.1 R_200 and R_200")
    #produce_single_simulation_graphs([set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types), 1.0), set_R_200_fraction(set_particle_types(sphere_specific_angular_momentum, *baryonic_low_mass_particle_types), 0.1)],
    #                                 ["R = $R_{200}$", "R = 0.1 $R_{200}$"], y_axis_label = "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Baryonic Specific Angular Momentum by particles within radius R at z=0",
    #                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #                                 filename = "j_baryon_inner_halo_and_halo_comparison.png")

    ## Dark Matter
    #print("Dark Matter - j - 0.1 R_200 / R_200")
    #produce_simulations_graph(lambda *args, **kwargs: sphere_specific_angular_momentum(particle_types = [ParticleType.dark_matter], R_200_fraction = 0.1, *args, **kwargs) / sphere_specific_angular_momentum(particle_types = [ParticleType.dark_matter], R_200_fraction = 1.0, *args, **kwargs),
    #                          "|$\\vec{j_{DM, 0.1 R_{200}}}$ / $\\vec{j_{DM, R_{200}}}$| (Arbitrary Units)", "Ratio of Dark Matter Specific Angular Momentum for particles within 0.1 $R_{200}$ to within $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = False, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1),
    #                          filename = "j_DM_inner_halo_to_halo_fraction.png")

    ## Baryons
    #print("Baryons - j - 0.1 R_200 / R_200")
    #produce_simulations_graph(lambda *args, **kwargs: sphere_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, R_200_fraction = 0.1, *args, **kwargs) / sphere_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, R_200_fraction = 1.0, *args, **kwargs),
    #                          "|$\\vec{j_{b, 0.1 R_{200}}}$ / $\\vec{j_{b, R_{200}}}$| (Arbitrary Units)", "Ratio of Baryonic Specific Angular Momentum for particles within 0.1 $R_{200}$ to within $R_{200}$ at z=0",
    #                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = False, invert_x = False, use_rolling_average = False,
    #                          xlim_overide = (0.1, 1),
    #                          filename = "j_baryon_inner_halo_to_halo_fraction.png")



    #for r_200_fraction in (1.0, 0.1):
    #    max_r_tag = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
    #    max_r_value = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
    #    for tag in constants.tags:
    #        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
    #            halo = assembily_history[sim][tag]["halo"]
    #            subhalo = assembily_history[sim][tag]["subhalo"]
    #            if halo is None:
    #                continue
    #            try:
    #                r = sphere_radius(halo, subhalo, tag, sim, [ParticleType.dark_matter], r_200_fraction)
    #            except LookupError:
    #                continue
    #            if r >= max_r_value[sim]:
    #                max_r_value[sim] = r
    #                max_r_tag[sim] = tag

    #    for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
    #        halo = assembily_history[sim][max_r_tag[sim]]["halo"]
    #        subhalo = assembily_history[sim][max_r_tag[sim]]["subhalo"]
    #        r_max_j = sphere_specific_angular_momentum(halo, subhalo, max_r_tag[sim], sim, [ParticleType.dark_matter], r_200_fraction)

    #        halo = assembily_history[sim][constants.tags[-1]]["halo"]
    #        subhalo = assembily_history[sim][constants.tags[-1]]["subhalo"]
    #        final_j = sphere_specific_angular_momentum(halo, subhalo, constants.tags[-1], sim, [ParticleType.dark_matter], r_200_fraction)

    #        ratio = final_j / r_max_j

    #        r_max_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", max_r_tag[sim], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][max_r_tag[sim]]["halo"] - 1]
    #        final_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", constants.tags[-1], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][constants.tags[-1]]["halo"] - 1]

    #        print(f"--|| INFO ||-- {sim} at R = {r_200_fraction} for Dark Matter:")
    #        print(f"                   j ratio     = {ratio}")
    #        print(f"                   r max tag   = {max_r_tag[sim]}")
    #        print(f"                   r max       = {max_r_value[sim]} kpc")
    #        print(f"                   r max j     = {r_max_j} kpc km / s")
    #        print(f"                   final j     = {final_j} kpc km / s")
    #        print(f"                   r max R_200 = {r_max_r_200} cm")
    #        print(f"                   final R_200 = {final_r_200} cm")
    #        print()



    for r_200_fraction in (1.0, 0.1):
        max_j_tag = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
        max_j_value = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
        for tag in constants.tags:
            for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
                halo = assembily_history[sim][tag]["halo"]
                subhalo = assembily_history[sim][tag]["subhalo"]
                if halo is None:
                    continue
                try:
                    j = sphere_specific_angular_momentum(halo, subhalo, tag, sim, [ParticleType.dark_matter], r_200_fraction)
                except LookupError:
                    continue
                if j >= max_j_value[sim]:
                    max_j_value[sim] = j
                    max_j_tag[sim] = tag

        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
            halo = assembily_history[sim][max_j_tag[sim]]["halo"]
            subhalo = assembily_history[sim][max_j_tag[sim]]["subhalo"]
            max_j = sphere_specific_angular_momentum(halo, subhalo, max_j_tag[sim], sim, [ParticleType.dark_matter], r_200_fraction)

            halo = assembily_history[sim][constants.tags[-1]]["halo"]
            subhalo = assembily_history[sim][constants.tags[-1]]["subhalo"]
            final_j = sphere_specific_angular_momentum(halo, subhalo, constants.tags[-1], sim, [ParticleType.dark_matter], r_200_fraction)

            ratio = final_j / max_j

            j_max_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", max_j_tag[sim], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][max_j_tag[sim]]["halo"] - 1]
            final_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", constants.tags[-1], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][constants.tags[-1]]["halo"] - 1]

            print(f"--|| INFO ||-- {sim} at R = {r_200_fraction} for Dark Matter:")
            print(f"                   j ratio     = {ratio}")
            print(f"                   j max tag   = {max_j_tag[sim]}")
            print(f"                   j max j     = {max_j} kpc km / s")
            print(f"                   final j     = {final_j} kpc km / s")
            print(f"                   j max R_200 = {j_max_r_200} cm")
            print(f"                   final R_200 = {final_r_200} cm")
            print()



    #for r_200_fraction in (1.0, 0.1):
    #    max_r_tag = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
    #    max_r_value = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
    #    for tag in constants.tags:
    #        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
    #            halo = assembily_history[sim][tag]["halo"]
    #            subhalo = assembily_history[sim][tag]["subhalo"]
    #            if halo is None:
    #                continue
    #            try:
    #                r = sphere_radius(halo, subhalo, tag, sim, [ParticleType.star], r_200_fraction)
    #            except LookupError:
    #                continue
    #            if r >= max_r_value[sim]:
    #                max_r_value[sim] = r
    #                max_r_tag[sim] = tag

    #    for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
    #        halo = assembily_history[sim][max_r_tag[sim]]["halo"]
    #        subhalo = assembily_history[sim][max_r_tag[sim]]["subhalo"]
    #        r_max_j = sphere_specific_angular_momentum(halo, subhalo, max_r_tag[sim], sim, [ParticleType.star], r_200_fraction)

    #        halo = assembily_history[sim][constants.tags[-1]]["halo"]
    #        subhalo = assembily_history[sim][constants.tags[-1]]["subhalo"]
    #        final_j = sphere_specific_angular_momentum(halo, subhalo, constants.tags[-1], sim, [ParticleType.star], r_200_fraction)

    #        ratio = final_j / r_max_j

    #        r_max_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", max_r_tag[sim], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][max_r_tag[sim]]["halo"] - 1]
    #        final_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", constants.tags[-1], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][constants.tags[-1]]["halo"] - 1]

    #        print(f"--|| INFO ||-- {sim} at R = {r_200_fraction} for Stars:")
    #        print(f"                   j ratio     = {ratio}")
    #        print(f"                   r max tag   = {max_r_tag[sim]}")
    #        print(f"                   r max       = {max_r_value[sim]} kpc")
    #        print(f"                   r max j     = {r_max_j} kpc km / s")
    #        print(f"                   final j     = {final_j} kpc km / s")
    #        print(f"                   r max R_200 = {r_max_r_200} cm")
    #        print(f"                   final R_200 = {final_r_200} cm")
    #        print()



    for r_200_fraction in (1.0, 0.1):
        max_j_tag = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
        max_j_value = { Simulations.Early:-1, Simulations.Organic:-1, Simulations.Late:-1 }
        for tag in constants.tags:
            for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
                halo = assembily_history[sim][tag]["halo"]
                subhalo = assembily_history[sim][tag]["subhalo"]
                if halo is None:
                    continue
                try:
                    j = sphere_specific_angular_momentum(halo, subhalo, tag, sim, [ParticleType.star], r_200_fraction)
                except LookupError:
                    continue
                if j >= max_j_value[sim]:
                    max_j_value[sim] = j
                    max_j_tag[sim] = tag

        for sim in (Simulations.Early, Simulations.Organic, Simulations.Late):
            halo = assembily_history[sim][max_j_tag[sim]]["halo"]
            subhalo = assembily_history[sim][max_j_tag[sim]]["subhalo"]
            max_j = sphere_specific_angular_momentum(halo, subhalo, max_j_tag[sim], sim, [ParticleType.star], r_200_fraction)

            halo = assembily_history[sim][constants.tags[-1]]["halo"]
            subhalo = assembily_history[sim][constants.tags[-1]]["subhalo"]
            final_j = sphere_specific_angular_momentum(halo, subhalo, constants.tags[-1], sim, [ParticleType.star], r_200_fraction)

            ratio = final_j / max_j

            j_max_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", max_j_tag[sim], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][max_j_tag[sim]]["halo"] - 1]
            final_r_200 = load_catalouge_field("Group_R_Crit200", "FOF", constants.tags[-1], sim, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[sim][constants.tags[-1]]["halo"] - 1]

            print(f"--|| INFO ||-- {sim} at R = {r_200_fraction} for Stars:")
            print(f"                   j ratio     = {ratio}")
            print(f"                   j max tag   = {max_j_tag[sim]}")
            print(f"                   j max j     = {max_j} kpc km / s")
            print(f"                   final j     = {final_j} kpc km / s")
            print(f"                   j max R_200 = {j_max_r_200} cm")
            print(f"                   final R_200 = {final_r_200} cm")
            print()
