import numpy as np
from matplotlib import pyplot as plt
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, barionic_particle_types, barionic_low_mass_particle_types, UnitSystem
from assembily_history import physical_centeral_mass_positions, assembily_history
from z0_particle_IDs import load_particle_IDs
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs, X_Axis_Value
from Physics import specific_angular_momentum, centre_of_mass
from angular_momentum_units import specific_angular_momentum_units
#TODO: try normalising by value of R_200?????

relitive_data_root = ".\\gm_for_mphys"

R_200_fraction = 1.0
#R_200_fraction = 0.1

def set_particle_types(func, *particle_types):
    def wrapper(*args, **kwargs):
        return func(particle_types = particle_types, *args, **kwargs)
    return wrapper

def r200_specific_angular_momentum(halo, subhalo, tag, simulation, particle_types = all_particle_types):
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
                
    #all_positions = np.empty(shape = (0, 3))
    #all_velocities = np.empty(shape = (0, 3))
    #all_masses = np.empty(shape = (0, ))
    #for particle_type in particle_IDs_by_type:
    #    box_particle_IDs, box_particle_locations = snapshot.particle_read(particle_type, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)
    #    id_filter = np.in1d(box_particle_IDs, relivant_particle_IDs)
    #    if sum(id_filter) == 0:
    #        continue# No particles with matching IDs present in this particle type
    #    #selection_radius = max(np.sqrt(((box_particle_locations[id_filter] - galaxy_centre)**2).sum(axis = 1))) + 0.001# Add a little bit to account for selection of smaler radi only
    #    selection_radius = max(np.sqrt(((box_particle_locations[id_filter] - selection_centre)**2).sum(axis = 1))) + 0.001# Add a little bit to account for selection of smaler radi only
    #    
    #    #particle_IDs, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "ParticleIDs", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    #    particle_IDs, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "ParticleIDs", selection_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    #    id_filter = np.in1d(particle_IDs, relivant_particle_IDs)
    #    particle_IDs = particle_IDs[id_filter]
    #    particle_locations_box_adjusted = particle_locations_box_adjusted[id_filter]
    #    #particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)[id_filter]
    #    particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", selection_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)[id_filter]
    #    if particle_type != ParticleType.dark_matter:
    #        #particle_masses = snapshot.particle_read_sphere(particle_type, "Mass", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)[id_filter]
    #        particle_masses = snapshot.particle_read_sphere(particle_type, "Mass", selection_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)[id_filter]
    #    else:
    #        particle_masses = np.full(particle_locations_box_adjusted.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        


        # Adjust for centre of mass of selected particles
        #initial_positions = particle_locations_box_adjusted
        #masses = particle_masses
        ##com = np.sum(initial_positions * masses[:, None], axis = 0) / sum(masses)
        #com = centre_of_mass(initial_positions, masses)
        #radius_adjustment = np.sqrt(np.sum(com**2))
        #
        #particle_IDs, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "ParticleIDs", galaxy_centre + com, selection_radius + radius_adjustment, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        #id_filter = np.in1d(particle_IDs, relivant_particle_IDs)
        #particle_IDs = particle_IDs[id_filter]
        #particle_locations_box_adjusted = particle_locations_box_adjusted[id_filter]
        #particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", galaxy_centre + com, selection_radius + radius_adjustment, UnitSystem.cgs, UnitSystem.cgs)[id_filter]
        #if particle_type != ParticleType.dark_matter:
        #    particle_masses = snapshot.particle_read_sphere(particle_type, "Mass", galaxy_centre + com, selection_radius + radius_adjustment, UnitSystem.cgs, UnitSystem.cgs)[id_filter]
        #else:
        #    particle_masses = np.full(particle_locations_box_adjusted.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        

        
    #    all_positions = np.append(all_positions, particle_locations_box_adjusted, axis = 0)
    #    all_velocities = np.append(all_velocities, particle_velocities, axis = 0)
    #    all_masses = np.append(all_masses, particle_masses)
#
    #return specific_angular_momentum_units(np.sqrt((specific_angular_momentum(all_positions, all_velocities, all_masses)**2).sum()))



# Dark Matter
produce_simulations_graph(set_particle_types(r200_specific_angular_momentum, ParticleType.dark_matter),
                          "|$\\vec{j_{DM}}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Dark Matter particles within {}$R_{{200}}$ at z=0".format("" if R_200_fraction == 1.0 else "{} ".format(R_200_fraction)),
                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
                          ylim_overide = (500, 20000),
                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$"))

# Barions (no Black Holes)
#produce_simulations_graph(set_particle_types(r200_specific_angular_momentum, *barionic_low_mass_particle_types),
#                          "|$\\vec{j_b}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Barion Particles within {}$R_{{200}}$ at z=0 (excluding Black Holes)".format("" if R_200_fraction == 1.0 else "{} ".format(R_200_fraction)),
#                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
#                          extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$"))

# Dark Matter vs. Barions
#produce_single_simulation_graphs([set_particle_types(r200_specific_angular_momentum, ParticleType.dark_matter), set_particle_types(r200_specific_angular_momentum, *barionic_low_mass_particle_types)],
#                                 ["Dark Matter", "Barionic Matter"], y_axis_label = "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", graph_title_partial = "Specific Angular Momentum by Matter Type (particles within {}$R_{{200}}$ at z=0)".format("" if R_200_fraction == 1.0 else "{} ".format(R_200_fraction)),
#                                 x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False)

# Barions / Dark Matter
#produce_simulations_graph(lambda *args, **kwargs: r200_specific_angular_momentum(particle_types = barionic_low_mass_particle_types, *args, **kwargs) / r200_specific_angular_momentum(particle_types = [ParticleType.dark_matter], *args, **kwargs),
#                          "|$\\vec{j_b}$ / $\\vec{j_b}$| (Arbitrary Units)", "Ratio of Barionic to Dark Matter Specific Angular Momentum for particles within {}$R_{{200}}$ at z=0 (excluding Black Holes)".format("" if R_200_fraction == 1.0 else "{} ".format(R_200_fraction)),
#                          x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False)
