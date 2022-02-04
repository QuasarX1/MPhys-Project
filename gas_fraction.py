from DataAccess import constants, load_catalouge_field, ParticleReadConversion_EagleSnapshot, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem

from data_over_time import produce_simulations_graph

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
    gas_mass = snapshot.particle_read_sphere(ParticleType.gas, "Mass", centre_point, r200, UnitSystem.cgs).sum()

    #stellar_mass = snapshot.particle_read_sphere(ParticleType.star, "Mass", centre_point, r200, UnitSystem.cgs).sum()
    #black_hole_mass = snapshot.particle_read_sphere(ParticleType.black_hole, "Mass", centre_point, r200, UnitSystem.cgs).sum()
    #dark_matter_mass = UnitSystem.convert_data(snapshot.header["MassTable"][ParticleType.dark_matter.value] * snapshot.header["NumPart_Total"][ParticleType.dark_matter.value], UnitSystem.h_less_comoving_GADGET, UnitSystem.cgs, cgs_conversion_factor = 10**10 / snapshot.hubble_paramiter)
    
    return (snapshot.header["Omega0"] / snapshot.header["OmegaBaryon"]) * gas_mass / m_200

produce_simulations_graph(func, "$f_{CGM}$ x $\Omega_0$ / $\Omega_b$", "Normalised Gas Mass Fraction", log_x = True)
