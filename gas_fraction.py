from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem

from data_over_time import produce_simulations_graph, X_Axis_Value

relitive_data_root = ".\\gm_for_mphys"

def func(halo, subhalo, tag, simulation):
    first_subhalo_index = load_catalouge_field("FirstSubhaloID", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root)[halo - 1]
    subhalo_index = first_subhalo_index + subhalo

    # Get the potential centre in data units
    centre_point = load_catalouge_field("CentreOfPotential", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root)[subhalo_index]

    # Get the R_200 from the group/halo entry in data units
    r200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root)[halo - 1]
    m_200 = load_catalouge_field("Group_M_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[halo - 1]

    # Get massesfor each particle type in R_200
    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulation, SimulationModels.RECAL, relitive_data_root)
    gas_masses = snapshot.particle_read_sphere(ParticleType.gas, "Mass", centre_point, r200, UnitSystem.cgs)

    # Gas Fraction
    #gas_mass = snapshot.particle_read_sphere(ParticleType.gas, "Mass", centre_point, r200, UnitSystem.cgs).sum()

    # Circumgalactic Gas Fraction
    star_formation_filter = snapshot.particle_read_sphere(ParticleType.gas, "StarFormationRate", centre_point, r200, UnitSystem.cgs) == 0.0
    gas_mass = snapshot.particle_read_sphere(ParticleType.gas, "Mass", centre_point, r200, UnitSystem.cgs)[star_formation_filter].sum()

    return (snapshot.header["Omega0"] / snapshot.header["OmegaBaryon"]) * gas_mass / m_200

produce_simulations_graph(func, "$f_{CGM}$ x $\Omega_0$ / $\Omega_b$", "Normalised Gas Mass Fraction",
                          x_axis = X_Axis_Value.time, log_x = False, invert_x = False,
                              xlim_overide = (0, constants.times[constants.tags[-1]]),
                              filename = "f_CGM_R200.png")
