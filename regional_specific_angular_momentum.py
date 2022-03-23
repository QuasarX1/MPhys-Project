import numpy as np
from matplotlib import pyplot as plt
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, barionic_particle_types, UnitSystem
from assembily_history import physical_centeral_mass_positions, assembily_history
from z0_particle_IDs import load_particle_IDs
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs, X_Axis_Value
from Physics import specific_angular_momentum
from angular_momentum_units import specific_angular_momentum_units
#TODO: try normalising by value of R_200?????

relitive_data_root = ".\\gm_for_mphys"

def set_particle_types(func, *particle_types):
    def wrapper(*args, **kwargs):
        return func(particle_types = particle_types, *args, **kwargs)
    return wrapper

def r200_specific_angular_momentum(halo, subhalo, tag, simulation, particle_types = all_particle_types):
    print("Doing {} {}".format(tag, Simulations.to_string(simulation)))
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    particle_IDs_by_type = load_particle_IDs(simulation, 1.0)
    relivant_particle_IDs = np.empty((sum([len(particle_IDs_by_type[particle_type]) for particle_type in particle_types]), ))
    next_index = 0
    for particle_type in particle_types:
        relivant_particle_IDs[next_index : next_index + len(particle_IDs_by_type[particle_type])] = particle_IDs_by_type[particle_type]
        next_index += len(particle_IDs_by_type[particle_type])

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)
                
    net_specific_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    for particle_type in particle_IDs_by_type:
        box_particle_IDs, box_particle_locations = snapshot.particle_read(particle_type, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)
        id_filter = np.in1d(box_particle_IDs, relivant_particle_IDs)
        if sum(id_filter) == 0:
            continue# No particles with matching IDs present in this particle type
        selection_radius = max(np.sqrt(((box_particle_locations[id_filter] - galaxy_centre)**2).sum(axis = 1))) + 0.001# Add a little bit to account for selection of smaler radi only
        
        particle_IDs, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "ParticleIDs", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        id_filter = np.in1d(particle_IDs, relivant_particle_IDs)
        particle_IDs = particle_IDs[id_filter]
        particle_locations_box_adjusted = particle_locations_box_adjusted[id_filter]
        particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)[id_filter]

        if particle_type != ParticleType.dark_matter:
            particle_masses = snapshot.particle_read_sphere(particle_type, "Mass", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)[id_filter]
        else:
            particle_masses = np.full(particle_locations_box_adjusted.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        
        #TODO: the return value might be a mass weighted average but surely this isnt!!!!!
        net_specific_angular_momentum += specific_angular_momentum(particle_locations_box_adjusted, particle_velocities, particle_masses)

    return specific_angular_momentum_units(np.sqrt((net_specific_angular_momentum**2).sum()))


#def alpha_exp_line(x_values, *args, **kwargs):
#    plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$")


#produce_simulations_graph(set_particle_types(r200_specific_angular_momentum, ParticleType.dark_matter),
#                          "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within $R_{200}$ at z=0",
#                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
#                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$"))

produce_simulations_graph(set_particle_types(r200_specific_angular_momentum, *barionic_particle_types),
                          "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Barion particles within $R_{200}$ at z=0",
                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**6, label = "j $\\propto$ $\\alpha^{3/2}$"))

produce_single_simulation_graphs([set_particle_types(r200_specific_angular_momentum, ParticleType.dark_matter), set_particle_types(r200_specific_angular_momentum, *barionic_particle_types)],
                                 ["Dark Matter", "Barionic Matter"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within $R_{200}$ at z=0)",
                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False)
