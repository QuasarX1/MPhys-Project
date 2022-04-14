from z0_ISM_passthrough_particle_IDs import load_particle_data
from z0_particle_IDs import load_particle_IDs as load_z0_R200_particle_IDs
from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, all_particle_types, baryonic_particle_types, baryonic_low_mass_particle_types, UnitSystem
import numpy as np
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs, X_Axis_Value, line_colours, line_styles
from Physics import specific_angular_momentum, centre_of_mass
from angular_momentum_units import specific_angular_momentum_units

relitive_data_root = ".\\gm_for_mphys"
SNe_feedback_temp_threshold = 10**7.5
AGN_feedback_temp_threshold = 10**8.5

def output_percentages(R_200_fraction = 1.0):
    for simulation in (Simulations.Early, Simulations.Organic, Simulations.Late):
        IDs, last_star_forming_a, max_temps, max_temp_a = load_particle_data(simulation, R_200_fraction)
        num_particles = len(IDs)
        print(f"--|| INFO ||-- Loaded {num_particles} particles from {simulation} within a radius of {R_200_fraction} * R_200")

        particles_experienced_feedback_filter = [max_temps >= SNe_feedback_temp_threshold]
        particles_experienced_SNe_feedback_filter = [(max_temps >= SNe_feedback_temp_threshold) & (max_temps <= AGN_feedback_temp_threshold)]
        particles_experienced_AGN_feedback_filter = [max_temps >= AGN_feedback_temp_threshold]
        print("--|| INFO ||-- Calculated filters")

        print("Feedback percentage = {}".format(np.sum(particles_experienced_feedback_filter) * 100 / num_particles))
        print("SNe only percentage = {}".format(np.sum(particles_experienced_SNe_feedback_filter) * 100 / num_particles))
        print("AGN only percentage = {}".format(np.sum(particles_experienced_AGN_feedback_filter) * 100 / num_particles))
        print("--|| INFO ||-- Done")
        print()
        print()
        print()

def set_R_200_fraction(func, fraction):
    def wrapper(*args, **kwargs):
        return func(R_200_fraction = fraction, *args, **kwargs)
    return wrapper

def sphere_specific_angular_momentum(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    print("--|| INFO ||-- Doing {} {}".format(tag, Simulations.to_string(simulation)))
    IDs, last_star_forming_a, max_temps, max_temp_a = load_particle_data(simulation, R_200_fraction)
    particleType = ParticleType.gas

    particles_experienced_feedback_filter = [max_temps >= SNe_feedback_temp_threshold]
    relivant_particle_IDs = IDs[particles_experienced_feedback_filter]

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    box_particle_IDs, box_particle_locations = snapshot.particle_read(particleType, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)

    id_filter = np.in1d(box_particle_IDs, relivant_particle_IDs)
    selected_box_particle_locations = box_particle_locations[id_filter]
    selected_box_particle_velocities = snapshot.particle_read(particleType, "Velocity", unit_system = UnitSystem.cgs)[id_filter]
    selected_box_particle_masses = snapshot.particle_read(particleType, "Mass", unit_system = UnitSystem.cgs)[id_filter]
        
    # Get R_200 set and use that for the centre
    R200_particle_IDs_by_type = load_z0_R200_particle_IDs(simulation)
    particle_types = all_particle_types
    R200_relivant_particle_IDs = np.empty((sum([len(R200_particle_IDs_by_type[particle_type]) for particle_type in particle_types]), ))
    next_index = 0
    for particle_type in particle_types:
        R200_relivant_particle_IDs[next_index : next_index + len(R200_particle_IDs_by_type[particle_type])] = R200_particle_IDs_by_type[particle_type]
        next_index += len(R200_particle_IDs_by_type[particle_type])
    box_all_positions = np.empty(shape = (0, 3))
    box_all_velocities = np.empty(shape = (0, 3))
    box_all_masses = np.empty(shape = (0, ))
    for particle_type in R200_particle_IDs_by_type:
        R200_box_particle_IDs, R200_box_particle_locations = snapshot.particle_read(particle_type, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)
        id_filter = np.in1d(R200_box_particle_IDs, R200_relivant_particle_IDs)
        if sum(id_filter) == 0:
            continue# No particles with matching IDs present in this particle type
        R200_box_particle_locations = R200_box_particle_locations[id_filter]
        R200_box_particle_velocities = snapshot.particle_read(particle_type, "Velocity", unit_system = UnitSystem.cgs)[id_filter]
        if particle_type != ParticleType.dark_matter:
            R200_box_particle_masses = snapshot.particle_read(particle_type, "Mass", unit_system = UnitSystem.cgs)[id_filter]
        else:
            R200_box_particle_masses = np.full(R200_box_particle_locations.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        box_all_positions = np.append(box_all_positions, R200_box_particle_locations, axis = 0)
        box_all_velocities = np.append(box_all_velocities, R200_box_particle_velocities, axis = 0)
        box_all_masses = np.append(box_all_masses, R200_box_particle_masses)
    selection_centre = centre_of_mass(box_all_positions, box_all_masses)
    
    all_positions = snapshot.convert_distance_values(
                        snapshot.centre_particles(
                            snapshot.convert_distance_values(selected_box_particle_locations, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET),
                            snapshot.convert_distance_values(selection_centre, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET)
                        ),
                        UnitSystem.h_less_comoving_GADGET, UnitSystem.cgs
                    )

    return specific_angular_momentum_units(np.sqrt((specific_angular_momentum(all_positions, selected_box_particle_velocities, selected_box_particle_masses)**2).sum()))

def sphere_median_radius(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    print("--|| INFO ||-- Doing {} {}".format(tag, Simulations.to_string(simulation)))
    IDs, last_star_forming_a, max_temps, max_temp_a = load_particle_data(simulation, R_200_fraction)
    particleType = ParticleType.gas

    particles_experienced_feedback_filter = [max_temps >= SNe_feedback_temp_threshold]
    relivant_particle_IDs = IDs[particles_experienced_feedback_filter]

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    box_particle_IDs, box_particle_locations = snapshot.particle_read(particleType, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)

    id_filter = np.in1d(box_particle_IDs, relivant_particle_IDs)
    selected_box_particle_locations = box_particle_locations[id_filter]
        
    # Get R_200 set and use that for the centre
    R200_particle_IDs_by_type = load_z0_R200_particle_IDs(simulation)
    particle_types = all_particle_types
    R200_relivant_particle_IDs = np.empty((sum([len(R200_particle_IDs_by_type[particle_type]) for particle_type in particle_types]), ))
    next_index = 0
    for particle_type in particle_types:
        R200_relivant_particle_IDs[next_index : next_index + len(R200_particle_IDs_by_type[particle_type])] = R200_particle_IDs_by_type[particle_type]
        next_index += len(R200_particle_IDs_by_type[particle_type])
    box_all_positions = np.empty(shape = (0, 3))
    box_all_masses = np.empty(shape = (0, ))
    for particle_type in R200_particle_IDs_by_type:
        R200_box_particle_IDs, R200_box_particle_locations = snapshot.particle_read(particle_type, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)
        id_filter = np.in1d(R200_box_particle_IDs, R200_relivant_particle_IDs)
        if sum(id_filter) == 0:
            continue# No particles with matching IDs present in this particle type
        R200_box_particle_locations = R200_box_particle_locations[id_filter]
        if particle_type != ParticleType.dark_matter:
            R200_box_particle_masses = snapshot.particle_read(particle_type, "Mass", unit_system = UnitSystem.cgs)[id_filter]
        else:
            R200_box_particle_masses = np.full(R200_box_particle_locations.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        box_all_positions = np.append(box_all_positions, R200_box_particle_locations, axis = 0)
        box_all_masses = np.append(box_all_masses, R200_box_particle_masses)
    selection_centre = centre_of_mass(box_all_positions, box_all_masses)
    
    all_positions = snapshot.convert_distance_values(
                        snapshot.centre_particles(
                            snapshot.convert_distance_values(selected_box_particle_locations, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET),
                            snapshot.convert_distance_values(selection_centre, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET)
                        ),
                        UnitSystem.h_less_comoving_GADGET, UnitSystem.cgs
                    )
    radi = np.sqrt(np.sum(all_positions**2, axis = 1))
    
    return np.median(radi) * 10**(3) / constants.SimulationConstants.get_constants()["CM_PER_MPC"]











def sphere_test_metric(halo, subhalo, tag, simulation, R_200_fraction = 1.0):
    print("--|| INFO ||-- Doing {} {}".format(tag, Simulations.to_string(simulation)))
    IDs, last_star_forming_a, max_temps, max_temp_a = load_particle_data(simulation, R_200_fraction)
    particleType = ParticleType.gas

    particles_experienced_feedback_filter = [max_temps >= SNe_feedback_temp_threshold]
    relivant_particle_IDs = IDs[particles_experienced_feedback_filter]

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    box_particle_IDs, box_particle_locations = snapshot.particle_read(particleType, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)

    id_filter = np.in1d(box_particle_IDs, relivant_particle_IDs)
    selected_box_particle_locations = box_particle_locations[id_filter]
    selected_box_particle_velocities = snapshot.particle_read(particleType, "Velocity", unit_system = UnitSystem.cgs)[id_filter]
    selected_box_particle_masses = snapshot.particle_read(particleType, "Mass", unit_system = UnitSystem.cgs)[id_filter]
        
    # Get R_200 set and use that for the centre
    R200_particle_IDs_by_type = load_z0_R200_particle_IDs(simulation)
    particle_types = all_particle_types
    R200_relivant_particle_IDs = np.empty((sum([len(R200_particle_IDs_by_type[particle_type]) for particle_type in particle_types]), ))
    next_index = 0
    for particle_type in particle_types:
        R200_relivant_particle_IDs[next_index : next_index + len(R200_particle_IDs_by_type[particle_type])] = R200_particle_IDs_by_type[particle_type]
        next_index += len(R200_particle_IDs_by_type[particle_type])
    box_all_positions = np.empty(shape = (0, 3))
    box_all_velocities = np.empty(shape = (0, 3))
    box_all_masses = np.empty(shape = (0, ))
    for particle_type in R200_particle_IDs_by_type:
        R200_box_particle_IDs, R200_box_particle_locations = snapshot.particle_read(particle_type, "ParticleIDs", unit_system = UnitSystem.cgs, return_coordinates = True)
        id_filter = np.in1d(R200_box_particle_IDs, R200_relivant_particle_IDs)
        if sum(id_filter) == 0:
            continue# No particles with matching IDs present in this particle type
        R200_box_particle_locations = R200_box_particle_locations[id_filter]
        R200_box_particle_velocities = snapshot.particle_read(particle_type, "Velocity", unit_system = UnitSystem.cgs)[id_filter]
        if particle_type != ParticleType.dark_matter:
            R200_box_particle_masses = snapshot.particle_read(particle_type, "Mass", unit_system = UnitSystem.cgs)[id_filter]
        else:
            R200_box_particle_masses = np.full(R200_box_particle_locations.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
        box_all_positions = np.append(box_all_positions, R200_box_particle_locations, axis = 0)
        box_all_velocities = np.append(box_all_velocities, R200_box_particle_velocities, axis = 0)
        box_all_masses = np.append(box_all_masses, R200_box_particle_masses)
    selection_centre = centre_of_mass(box_all_positions, box_all_masses)
    
    all_positions = snapshot.convert_distance_values(
                        snapshot.centre_particles(
                            snapshot.convert_distance_values(selected_box_particle_locations, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET),
                            snapshot.convert_distance_values(selection_centre, UnitSystem.cgs, UnitSystem.h_less_comoving_GADGET)
                        ),
                        UnitSystem.h_less_comoving_GADGET, UnitSystem.cgs
                    )

    radi = np.sqrt(np.sum(all_positions**2, axis = 1))
    
    median_radius = np.median(radi) * 10**(3) / constants.SimulationConstants.get_constants()["CM_PER_MPC"]

    return specific_angular_momentum_units(np.sqrt((specific_angular_momentum(all_positions, selected_box_particle_velocities, selected_box_particle_masses)**2).sum())) / median_radius









if __name__ == "__main__":
    #R_200_fraction = 0.1
    #R_200_fraction = 1.0
    #R_200_fraction = 3.0
    #R_200_fraction = None
    #output_percentages(R_200_fraction)

    # Median R
    print("Gas - Median r - whole box")
    produce_simulations_graph(set_R_200_fraction(sphere_median_radius, None),
        "Median Population Radius ($kPc$)", "Median Radius of the z=0 population of Gas Particles that have passed through the ISM and Experienced Feedback",
        x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
        xlim_overide = (0.1, 1), ylim_overide = (100, 3000),
        filename = "r_median_feedback_ISM_gas.png")
    print("--|| INF0 ||-- Done")
    
    # Median R at 3.0 * R_200
    print("Gas - Median r - 3.0 R_200")
    produce_simulations_graph(set_R_200_fraction(sphere_median_radius, 3.0),
        "Median Population Radius ($kPc$)", "Median Radius of the z=0 population of Gas Particles that have passed through the ISM and Experienced Feedback (R = 3.0 * R_200)",
        x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
        xlim_overide = (0.1, 1), ylim_overide = (100, 3000),
        filename = "r_median_feedback_ISM_gas_3.0_times_R_200.png")
    print("--|| INF0 ||-- Done")

    ## j
    #print("Gas - j - whole box")
    #produce_simulations_graph(set_R_200_fraction(sphere_specific_angular_momentum, None),
    #    "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Gas Particles that have passed through the ISM and Experienced Feedback",
    #    x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #    xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #    #extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]),
    #    filename = "j_feedback_ISM_gas.png")
    #print("--|| INF0 ||-- Done")

    ## j at 3.0 * R_200
    #print("Gas - j - 3.0 R_200")
    #produce_simulations_graph(set_R_200_fraction(sphere_specific_angular_momentum, 3.0),
    #    "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", "Specific Angular Momentum of Gas Particles that have passed through the ISM and Experienced Feedback (R = 3.0 * R_200)",
    #    x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #    xlim_overide = (0.1, 1), ylim_overide = (100, 20000),
    #    #extra_plotting_func = lambda x_values, *args, **kwargs: plt.plot(x_values[1][:5], kwargs["expansion_factor_values"][1][:5]**(3/2) * 10**5, label = "j $\\propto$ $\\alpha^{3/2}$", color = line_colours[-1], linestyle = line_styles[-1]),
    #    filename = "j_feedback_ISM_gas_3.0_times_R_200.png")
    #print("--|| INF0 ||-- Done")






    ## j/median(r)
    #print("Gas - j / median(r) - whole box")
    #produce_simulations_graph(set_R_200_fraction(sphere_test_metric, None),
    #    "|$\\vec{j}$| / median(r) ($km$ $s^{-1}$)", "Specific Angular Momentum of Gas Particles that have passed through the ISM and Experienced Feedback",
    #    x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #    xlim_overide = (0.1, 1), ylim_overide = (5, 60),
    #    filename = "TEST_j_over_r_feedback_ISM_gas.png")
    #print("--|| INF0 ||-- Done")

    ## j/median(r) at 3.0 * R_200
    #print("Gas - j / median(r) - 3.0 R_200")
    #produce_simulations_graph(set_R_200_fraction(sphere_test_metric, 3.0),
    #    "|$\\vec{j}$| / median(r) ($km$ $s^{-1}$)", "Specific Angular Momentum of Gas Particles that have passed through the ISM and Experienced Feedback (R = 3.0 * R_200)",
    #    x_axis = X_Axis_Value.expansion_factor, log_x = True, log_y = True, invert_x = False, use_rolling_average = False,
    #    xlim_overide = (0.1, 1), ylim_overide = (5, 60),
    #    filename = "TEST_j_over_r_feedback_ISM_gas_3.0_times_R_200.png")
    #print("--|| INF0 ||-- Done")

    
