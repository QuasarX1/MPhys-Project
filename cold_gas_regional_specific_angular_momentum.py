import numpy as np
from matplotlib import pyplot as plt
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, baryonic_particle_types, baryonic_low_mass_particle_types, UnitSystem
from assembily_history import physical_centeral_mass_positions, assembily_history
from cold_gas_z0_particle_IDs import load_particle_IDs
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

def r200_specific_angular_momentum(halo, subhalo, tag, simulation, particle_types = all_particle_types, R_200_fraction = 1.0):
    print("Doing {} {}".format(tag, Simulations.to_string(simulation)))
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")


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



# 1.0 R_200 ----------------------------------------------------------------

# Baryons (no Black Holes)
produce_simulations_graph(set_particle_types(r200_specific_angular_momentum, *baryonic_low_mass_particle_types),
                          "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Baryon Particles within $R_{200}$ at z=0 (excluding Black Holes)",
                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                          xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]))

# Dark Matter vs. Baryons
#produce_single_simulation_graphs([set_particle_types(r200_specific_angular_momentum, ParticleType.dark_matter), set_particle_types(r200_specific_angular_momentum, *baryonic_low_mass_particle_types)],
#                                 ["Dark Matter", "Baryonic Matter"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within $R_{200}$ at z=0)",
#                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
#                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000))

# Baryons / Dark Matter
#produce_simulations_graph(lambda *args, **kwargs: r200_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, *args, **kwargs) / r200_specific_angular_momentum(particle_types = [ParticleType.dark_matter], *args, **kwargs),
#                          "|$\\vec{j_b}$ / $\\vec{j_DM}$| (Arbitrary Units)", "Ratio of Baryonic to Dark Matter Specific Angular Momentum for particles within $R_{200}$ at z=0 (excluding Black Holes)",
#                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = False, invert_x = False, use_rolling_average = False,
#                          xlim_overide = (0.1, 1))



# 0.1 R_200 ----------------------------------------------------------------

# Baryons (no Black Holes)
produce_simulations_graph(set_R_200_fraction(set_particle_types(r200_specific_angular_momentum, *baryonic_low_mass_particle_types), 0.1),
                          "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Baryon Particles within 0.1 $R_{200}$ at z=0 (excluding Black Holes)",
                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                          xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]))

# Dark Matter vs. Baryons
#produce_single_simulation_graphs([set_R_200_fraction(set_particle_types(r200_specific_angular_momentum, ParticleType.dark_matter), 0.1), set_R_200_fraction(set_particle_types(r200_specific_angular_momentum, *baryonic_low_mass_particle_types), 0.1)],
#                                 ["Dark Matter", "Baryonic Matter"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within 0.1 $R_{200}$ at z=0)",
#                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
#                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000))

# Baryons / Dark Matter
#produce_simulations_graph(lambda *args, **kwargs: r200_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, R_200_fraction = 0.1, *args, **kwargs) / r200_specific_angular_momentum(particle_types = [ParticleType.dark_matter], R_200_fraction = 0.1, *args, **kwargs),
#                          "|$\\vec{j_b}$ / $\\vec{j_DM}$| (Arbitrary Units)", "Ratio of Baryonic to Dark Matter Specific Angular Momentum for particles within 0.1 $R_{200}$ at z=0 (excluding Black Holes)",
#                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = False, invert_x = False, use_rolling_average = False,
#                          xlim_overide = (0.1, 1))



# 0.1 R_200 / 1.0 R_200 ----------------------------------------------------

# Baryons
produce_simulations_graph(lambda *args, **kwargs: r200_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, R_200_fraction = 0.1, *args, **kwargs) / r200_specific_angular_momentum(particle_types = baryonic_low_mass_particle_types, R_200_fraction = 1.0, *args, **kwargs),
                          "|$\\vec{j_{b, 0.1 R_{200}}}$ / $\\vec{j_{b, R_{200}}}$| (Arbitrary Units)", "Ratio of Baryonic Specific Angular Momentum for particles within 0.1 $R_{200}$ to within $R_{200}$ at z=0",
                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = False, invert_x = False, use_rolling_average = False,
                          xlim_overide = (0.1, 1))

# Dark Matter Comparison
produce_single_simulation_graphs([set_R_200_fraction(set_particle_types(r200_specific_angular_momentum, *baryonic_low_mass_particle_types), 1.0), set_R_200_fraction(set_particle_types(r200_specific_angular_momentum, *baryonic_low_mass_particle_types), 0.1)],
                                 ["R = $R_{200}$", "R = 0.1 $R_{200}$"], y_axis_label = "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Baryonic Specific Angular Momentum by particles within radius R at z=0",
                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                                 xlim_overide = (0.1, 1), ylim_overide = (100, 20000))
