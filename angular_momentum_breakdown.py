import numpy as np

from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem

from assembily_history import physical_centeral_mass_positions, assembily_history
from data_over_time import produce_simulations_graph, produce_single_simulation_graphs

relitive_data_root = ".\\gm_for_mphys"

def set_particle_type(func, particle: ParticleType):
    def wrapper(*args, **kwargs):
        return func(particle, *args, **kwargs)
    return wrapper

def barionic_angular_momenta_data(halo, subhalo, tag, simulation):
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulation][tag]["halo"] - 1]
    #selection_radius = r_200
    selection_radius = 9.257*(10**22)# 30 KPc

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    # Calculate the net angular momentum (axis vector)
    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
    net_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    angular_momenta = {particle_type:np.array([0, 0, 0], dtype = np.float64) for particle_type in particle_types}
    for particle_type in particle_types:
        _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_type, "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_velocities = snapshot.particle_read_sphere(particle_type, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        particle_velocities[:, 0] -= particle_velocities[:, 0].mean()
        particle_velocities[:, 1] -= particle_velocities[:, 1].mean()
        particle_velocities[:, 2] -= particle_velocities[:, 2].mean()
        particle_masses = snapshot.particle_read_sphere(particle_type, "Mass", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)

        angular_momenta[particle_type] = np.sum(np.cross(particle_locations_box_adjusted, particle_velocities) * particle_masses[:, None], axis = 0)
        net_angular_momentum += angular_momenta[particle_type]

    return angular_momenta, net_angular_momentum

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

def gas_star_angular_momentum_angle(halo, subhalo, tag, simulation):
    angular_momenta, _ = barionic_angular_momenta_data(halo, subhalo, tag, simulation)

    return np.arccos(angular_momenta[ParticleType.gas].dot(angular_momenta[ParticleType.star]) / (np.linalg.norm(angular_momenta[ParticleType.gas], 2) * np.linalg.norm(angular_momenta[ParticleType.star], 2))) * 180 / np.pi

def DM_angular_momenta(selection_radius, halo, subhalo, tag, simulation):
    galaxy_centre = physical_centeral_mass_positions[simulation][tag]#TODO: check whether cgs or physical code is uncommented!!!
    if galaxy_centre is None:
        raise LookupError("No data avalible.")

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)

    _, particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.dark_matter, "Coordinates", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    particle_velocities = snapshot.particle_read_sphere(ParticleType.dark_matter, "Velocity", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    particle_velocities[:, 0] -= particle_velocities[:, 0].mean()
    particle_velocities[:, 1] -= particle_velocities[:, 1].mean()
    particle_velocities[:, 2] -= particle_velocities[:, 2].mean()
    #particle_masses = snapshot.particle_read_sphere(ParticleType.dark_matter, "Mass", galaxy_centre, selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    particle_masses = np.full(particle_locations_box_adjusted.shape[0], snapshot.header["MassTable"][ParticleType.dark_matter.value] * 10**10 / snapshot.hubble_paramiter)
    #particle_masses = UnitSystem.convert_data(snapshot.header["MassTable"][ParticleType.dark_matter.value] * snapshot.header["NumPart_Total"][ParticleType.dark_matter.value], UnitSystem.h_less_comoving_GADGET, UnitSystem.cgs, cgs_conversion_factor = 10**10 / snapshot.hubble_paramiter)

    angular_momenta = np.cross(particle_locations_box_adjusted, particle_velocities) * particle_masses[:, None]

    # Calculate the net angular momentum
    net_angular_momentum = np.sum(angular_momenta, axis = 0)

    return angular_momenta, net_angular_momentum

def DM_halo_angular_momenta(halo, subhalo, tag, simulation):
    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulation][tag]["halo"] - 1]
    return UnitSystem.convert_mass_to_solar(np.sqrt((DM_angular_momenta(r_200, halo, subhalo, tag, simulation)[1]**2).sum())) * 1000000 * constants.SimulationConstants.get_constants()["SEC_PER_MEGAYEAR"] / constants.SimulationConstants.get_constants()["CM_PER_MPC"]**2

def DM_30kpc_angular_momenta(halo, subhalo, tag, simulation):
    r = 9.257*(10**22)# 30 KPc
    return UnitSystem.convert_mass_to_solar(np.sqrt((DM_angular_momenta(r, halo, subhalo, tag, simulation)[1]**2).sum())) * 1000000 * constants.SimulationConstants.get_constants()["SEC_PER_MEGAYEAR"] / constants.SimulationConstants.get_constants()["CM_PER_MPC"]**2



produce_simulations_graph(set_particle_type(angular_momentum_fraction, ParticleType.gas), "Fraction of Total Galactic $\\vec{L}$", "Gaseous $\\vec{L}$ Fraction", log_x = True)
produce_simulations_graph(set_particle_type(angular_momentum_fraction, ParticleType.star), "Fraction of Total Galactic $\\vec{L}$", "Stellar $\\vec{L}$ Fraction", log_x = True)
produce_single_simulation_graphs([set_particle_type(angular_momentum_fraction, ParticleType.gas), set_particle_type(angular_momentum_fraction, ParticleType.star)], ["Gas", "Star"], "Fraction of Total Galactic $\\vec{L}$", "$\\vec{L}$ Fraction", log_x = True)
produce_simulations_graph(gas_star_angular_momentum_angle, "$\\theta$ (Degrees)", "Gaseous-Stellar $\\vec{L}$ Angular Seperation", log_x = True)

produce_simulations_graph(DM_halo_angular_momenta, "|$\\vec{L}$| ($M_{sun}$ $Pc$ $Myrs^{-1}$)", "Dark Matter ($R_{200}$) Angular Momentum", log_x = True)
produce_simulations_graph(DM_30kpc_angular_momenta, "|$\\vec{L}$| ($M_{sun}$ $Pc$ $Myrs^{-1}$)", "Dark Matter (30 KPc) Angular Momentum", log_x = True)
