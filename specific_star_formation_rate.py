from DataAccess import constants, load_catalouge_field, Simulations, SimulationModels, ParticleType, UnitSystem

from data_over_time import produce_simulations_graph, X_Axis_Value

relitive_data_root = ".\\gm_for_mphys"

def func(halo, subhalo, tag, simulation):
    first_subhalo_index = load_catalouge_field("FirstSubhaloID", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root)[halo - 1]
    subhalo_index = first_subhalo_index + subhalo

    stellar_mass = UnitSystem.convert_mass_to_solar(load_catalouge_field("MassType", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[subhalo_index][ParticleType.star.value])
    
    star_formation_rate = UnitSystem.convert_mass_to_solar(load_catalouge_field("StarFormationRate", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[subhalo_index])
    # Per second to per year convertion
    star_formation_rate *= constants.SimulationConstants.get_constants()["SEC_PER_YEAR"]

    return star_formation_rate / stellar_mass

produce_simulations_graph(func, "sSFR ($yr^{-1}$)", "Specific Star Formation Rate",
                          x_axis = X_Axis_Value.time, log_x = False, log_y = True, invert_x = False,
                          xlim_overide = (0, constants.times[constants.tags[-1]]),
                          filename = "sSFR_30kpc.png")
