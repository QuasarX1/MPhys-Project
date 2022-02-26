import numpy as np

from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem
from Physics import angular_momentum, specific_angular_momentum

from assembily_history import physical_centeral_mass_positions, assembily_history
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs

relitive_data_root = ".\\gm_for_mphys"

def set_particle_type(func, particle: ParticleType):
    def wrapper(*args, **kwargs):
        return func(particle, *args, **kwargs)
    return wrapper

def set_selection_radius(func, r200 = False):
    def wrapper(*args, **kwargs):
        if not r200:
            return func(*args, **kwargs)
        else:
            return func(use_r200 = True, *args, **kwargs)
    return wrapper


def angular_momentum_units(value):
    """
    g cm^2 s^-1 -> M_sun kPc km s^-1
    """
    return UnitSystem.convert_mass_to_solar(value) * 10**(-1) / constants.SimulationConstants.get_constants()["CM_PER_MPC"]


def specific_angular_momentum_units(value):
    """
    cm^2 s^-1 -> kPc km s^-1
    """
    return value * 10**(-1) / constants.SimulationConstants.get_constants()["CM_PER_MPC"]


def barionic_angular_momenta_data(halo, subhalo, tag, simulation, use_r200 = False):
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulation][tag]["halo"] - 1]
    if use_r200:
        selection_radius = r_200
    else:
        selection_radius = 9.257*(10**22)# 30 KPc

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
    net_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    angular_momenta = {particle_type:np.array([0, 0, 0], dtype = np.float64) for particle_type in particle_types}
    for particle_type in particle_types:
        # Read data
        _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        particle_masses = snapshot.particle_read_sphere(particle_type, "Mass", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        
        angular_momenta[particle_type] = angular_momentum(particle_locations_box_adjusted, particle_velocities, particle_masses)
        net_angular_momentum += angular_momenta[particle_type]

    return angular_momenta, net_angular_momentum


def barionic_specific_angular_momenta_data(halo, subhalo, tag, simulation, use_r200 = False):
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulation][tag]["halo"] - 1]
    if use_r200:
        selection_radius = r_200
    else:
        selection_radius = 9.257*(10**22)# 30 KPc

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
    net_specific_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    specific_angular_momenta = {particle_type:np.array([0, 0, 0], dtype = np.float64) for particle_type in particle_types}
    for particle_type in particle_types:
        # Read data
        _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        
        specific_angular_momenta[particle_type] = specific_angular_momentum(particle_locations_box_adjusted, particle_velocities)
        net_specific_angular_momentum += specific_angular_momenta[particle_type]

    return specific_angular_momenta, net_specific_angular_momentum


def angular_momentum_fraction(target_particle_type, halo, subhalo, tag, simulation):
    angular_momenta, net_angular_momentum = barionic_angular_momenta_data(halo, subhalo, tag, simulation)
        
    total_angular_momentum_length = np.linalg.norm(net_angular_momentum, 2)

    if isinstance(target_particle_type, (list, tuple, np.ndarray)):
        target_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
        for particle_type in target_particle_type:
            target_angular_momentum += angular_momenta[particle_type]

    else:
        target_angular_momentum = angular_momenta[target_particle_type]

    return net_angular_momentum.dot(target_angular_momentum) / (total_angular_momentum_length**2)


def specific_angular_momentum_fraction(target_particle_type, halo, subhalo, tag, simulation):
    specific_angular_momenta, net_specific_angular_momentum = barionic_specific_angular_momenta_data(halo, subhalo, tag, simulation)
        
    total_specific_angular_momentum_length = np.linalg.norm(net_specific_angular_momentum, 2)

    if isinstance(target_particle_type, (list, tuple, np.ndarray)):
        target_specific_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
        for particle_type in target_particle_type:
            target_specific_angular_momentum += specific_angular_momenta[particle_type]

    else:
        target_specific_angular_momentum = specific_angular_momenta[target_particle_type]

    return net_specific_angular_momentum.dot(target_specific_angular_momentum) / (total_specific_angular_momentum_length**2)


def gas_star_angular_momentum_angle(halo, subhalo, tag, simulation):
    angular_momenta, _ = barionic_angular_momenta_data(halo, subhalo, tag, simulation)

    return np.arccos(angular_momenta[ParticleType.gas].dot(angular_momenta[ParticleType.star]) / (np.linalg.norm(angular_momenta[ParticleType.gas], 2) * np.linalg.norm(angular_momenta[ParticleType.star], 2))) * 180 / np.pi


def DM_angular_momentum(selection_radius, halo, subhalo, tag, simulation):
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    _, particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.dark_matter, "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    particle_velocities = snapshot.particle_read_sphere(ParticleType.dark_matter, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    particle_masses = np.full(particle_locations_box_adjusted.shape[0], UnitSystem.convert_mass_from_solar(snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter))
    
    # Calculate the net angular momentum
    return angular_momentum(particle_locations_box_adjusted, particle_velocities, particle_masses)


def DM_halo_angular_momentum(halo, subhalo, tag, simulation):
    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulation][tag]["halo"] - 1]
    return angular_momentum_units(np.sqrt((DM_angular_momentum(r_200, halo, subhalo, tag, simulation)**2).sum()))


def DM_30kpc_angular_momentum(halo, subhalo, tag, simulation):
    r = 9.257*(10**22)# 30 KPc
    return angular_momentum_units(np.sqrt((DM_angular_momentum(r, halo, subhalo, tag, simulation)**2).sum()))


def DM_specific_angular_momentum(selection_radius, halo, subhalo, tag, simulation):
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    _, particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.dark_matter, "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    particle_velocities = snapshot.particle_read_sphere(ParticleType.dark_matter, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    
    # Calculate the net angular momentum
    return specific_angular_momentum(particle_locations_box_adjusted, particle_velocities)


def DM_halo_specific_angular_momentum(halo, subhalo, tag, simulation):
    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulation][tag]["halo"] - 1]
    return specific_angular_momentum_units(np.sqrt((DM_specific_angular_momentum(r_200, halo, subhalo, tag, simulation)**2).sum()))


def DM_30kpc_specific_angular_momentum(halo, subhalo, tag, simulation):
    r = 9.257*(10**22)# 30 KPc
    return specific_angular_momentum_units(np.sqrt((DM_specific_angular_momentum(r, halo, subhalo, tag, simulation)**2).sum()))


def particle_specific_angular_momentum(target_particle_type, halo, subhalo, tag, simulation, use_r200 = False):
    return np.linalg.norm(barionic_specific_angular_momenta_data(halo, subhalo, tag, simulation, use_r200)[0][target_particle_type])



#produce_simulations_graph(set_particle_type(angular_momentum_fraction, ParticleType.gas), "Fraction of Total Galactic $\\vec{L}$", "Gaseous $\\vec{L}$ Fraction", log_x = True)
#produce_simulations_graph(set_particle_type(angular_momentum_fraction, ParticleType.star), "Fraction of Total Galactic $\\vec{L}$", "Stellar $\\vec{L}$ Fraction", log_x = True)
#produce_single_simulation_graphs([set_particle_type(angular_momentum_fraction, ParticleType.gas), set_particle_type(angular_momentum_fraction, ParticleType.star)], ["Gas", "Star"], "Fraction of Total Galactic $\\vec{L}$", "$\\vec{L}$ Fraction", log_x = True)

#produce_simulations_graph(gas_star_angular_momentum_angle, "$\\theta$ (Degrees)", "Gaseous-Stellar $\\vec{L}$ Angular Seperation", log_x = True)

#produce_simulations_graph(DM_halo_angular_momentum, "|$\\vec{L}$| ($M_{sun}$ $kPc$ $km$ $s^{-1}$)", "Dark Matter ($R_{200}$) Angular Momentum", log_x = True)
#produce_simulations_graph(DM_30kpc_angular_momentum, "|$\\vec{L}$| ($M_{sun}$ $kPc$ $km$ $s^{-1}$)", "Dark Matter (30 KPc) Angular Momentum", log_x = True)


produce_simulations_graph(set_selection_radius(set_particle_type(particle_specific_angular_momentum, ParticleType.star), r200 = True), "Galactic |$\\vec{j_*}$| ($kPc$ $km$ $s^{-1}$)", "Stellar Specific Angular Momentum inside $R_{200}$", log_x = True)
produce_simulations_graph(set_selection_radius(set_particle_type(particle_specific_angular_momentum, ParticleType.gas), r200 = True), "Galactic |$\\vec{j_*}$| ($kPc$ $km$ $s^{-1}$)", "Gaseous Specific Angular Momentum inside $R_{200}$", log_x = True)
produce_simulations_graph(set_selection_radius(set_particle_type(particle_specific_angular_momentum, ParticleType.black_hole), r200 = True), "Galactic |$\\vec{j_*}$| ($kPc$ $km$ $s^{-1}$)", "Black Hole Specific Angular Momentum inside $R_{200}$", log_x = True)

#produce_simulations_graph(set_particle_type(specific_angular_momentum_fraction, ParticleType.gas), "Fraction of Total Galactic $\\vec{j}$", "Gaseous $\\vec{j}$ Fraction", log_x = True)
#produce_simulations_graph(set_particle_type(specific_angular_momentum_fraction, ParticleType.star), "Fraction of Total Galactic $\\vec{j}$", "Stellar $\\vec{j}$ Fraction", log_x = True)
#produce_single_simulation_graphs([set_particle_type(specific_angular_momentum_fraction, ParticleType.gas), set_particle_type(specific_angular_momentum_fraction, ParticleType.star)], ["Gas", "Star"], "Fraction of Total Galactic $\\vec{j}$", "$\\vec{j}$ Fraction", log_x = True)

produce_simulations_graph(DM_halo_specific_angular_momentum, "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", "Dark Matter ($R_{200}$) Specific Angular Momentum", log_x = True)
produce_simulations_graph(DM_30kpc_specific_angular_momentum, "|$\\vec{j}$| ($kPc$ $km$ $s^{-1}$)", "Dark Matter (30 KPc) Specific Angular Momentum", log_x = True)
