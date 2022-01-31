import numpy as np
from matplotlib import pyplot as plt

from DataAccess import constants, load_catalouge_field, particles_by_group, Simulations, SimulationModels, ParticleType, UnitSystem
from assembily_history import assembily_history

relitive_data_root = ".\\gm_for_mphys"

def total_mass_M200(halo, subhalo, tag, simulation):
    return UnitSystem.convert_mass_to_solar(load_catalouge_field("Group_M_Crit200", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[halo - 1])

def _particle_type_mass(particle_type, halo, subhalo, tag, simulation):
    first_subhalo_index = load_catalouge_field("FirstSubhaloID", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root)[halo - 1]
    subhalo_index = first_subhalo_index + subhalo
    return UnitSystem.convert_mass_to_solar(load_catalouge_field("ApertureMeasurements/Mass/030kpc", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[subhalo_index][particle_type.value])

def gas_mass(*args, **kwargs):
    return _particle_type_mass(ParticleType.gas, *args, **kwargs)

def stellar_mass(*args, **kwargs):
    return _particle_type_mass(ParticleType.star, *args, **kwargs)

def black_hole_mass(*args, **kwargs):
    return _particle_type_mass(ParticleType.black_hole, *args, **kwargs)

def dark_matter_mass(*args, **kwargs):
    return _particle_type_mass(ParticleType.dark_matter, *args, **kwargs)

def baryon_mass(halo, subhalo, tag, simulation):
    first_subhalo_index = load_catalouge_field("FirstSubhaloID", "FOF", tag, simulation, SimulationModels.RECAL, relitive_data_root)[halo - 1]
    subhalo_index = first_subhalo_index + subhalo
    masses = load_catalouge_field("ApertureMeasurements/Mass/030kpc", "Subhalo", tag, simulation, SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[subhalo_index]
    return UnitSystem.convert_mass_to_solar(masses[ParticleType.gas.value] + masses[ParticleType.star.value] + masses[ParticleType.black_hole.value])




def produce_simulations_graph(func, y_axis_label, graph_title_partial, log_x = False, log_y = False):
    mass_values = []
    snapshot_numbers = []
    for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
        mass_values.append(np.empty(shape = (len(constants.tags), )))
        snapshot_numbers.append(np.array([float("{}.{}".format(tag[5:8], tag[9:12])) for tag in constants.tags]))
        for j, tag in enumerate(constants.tags):
            try:
                halo = assembily_history[simulation][tag]["halo"]
                subhalo = assembily_history[simulation][tag]["subhalo"]
                if halo is None:
                    raise LookupError("No history avalible.")

                mass_values[i][j] = func(halo, subhalo, tag, simulation)

            except LookupError:
                mass_values[i][j] = None
            
        array_filter = True ^ np.isnan(mass_values[i])
        mass_values[i] = mass_values[i][array_filter]
        snapshot_numbers[i] = snapshot_numbers[i][array_filter]

        
    plt.plot(snapshot_numbers[0], mass_values[0], label = "Early")
    plt.plot(snapshot_numbers[1], mass_values[1], label = "Organic")
    plt.plot(snapshot_numbers[2], mass_values[2], label = "Late")
    if log_x and log_y:
        plt.loglog()
    elif log_x:
        plt.semilogx()
    elif log_y:
        plt.semilogy()
    plt.xlim(plt.xlim()[1], plt.xlim()[0])
    plt.xlabel("Redshift")
    plt.ylabel(y_axis_label)
    plt.title(f"{graph_title_partial} for each Snapshot")
    plt.legend()
    plt.show()

def produce_single_simulation_graphs(funcs, quantity_labels, y_axis_label, graph_title_partial, log_x = False, log_y = False):
    for simulation in (Simulations.Early, Simulations.Organic, Simulations.Late):
        quantity_values = []
        snapshot_numbers = []
        for i in range(len(funcs)):
            quantity_values.append(np.empty(shape = (len(constants.tags), )))
            snapshot_numbers.append(np.array([float("{}.{}".format(tag[5:8], tag[9:12])) for tag in constants.tags]))
        
        for j, tag in enumerate(constants.tags):
            try:
                halo = assembily_history[simulation][tag]["halo"]
                subhalo = assembily_history[simulation][tag]["subhalo"]
                if halo is None:
                    raise LookupError("No history avalible.")

                for i, func in enumerate(funcs):
                    try:
                        quantity_values[i][j] = func(halo, subhalo, tag, simulation)
                    except LookupError:
                        quantity_values[i][j] = None
                
            except LookupError:
                for i in range(len(funcs)):
                    quantity_values[i][j] = None
            
        for i in range(len(funcs)):
            array_filter = True ^ np.isnan(quantity_values[i])
            quantity_values[i] = quantity_values[i][array_filter]
            snapshot_numbers[i] = snapshot_numbers[i][array_filter]
        
            plt.plot(snapshot_numbers[i], quantity_values[i], label = quantity_labels[i])

        if log_x and log_y:
            plt.loglog()
        elif log_x:
            plt.semilogx()
        elif log_y:
            plt.semilogy()
        plt.xlim(plt.xlim()[1], plt.xlim()[0])
        plt.xlabel("Redshift")
        plt.ylabel(y_axis_label)
        simulation_name = Simulations.to_string(simulation)
        plt.title(f"{simulation_name} {graph_title_partial} for each Snapshot")
        plt.legend()
        plt.show()



if __name__ == "__main__":
    #produce_simulations_graph(total_mass_M200, y_axis_label = "$M_{200}$ ($M_{sun}$)", graph_title_partial = "Group_M_Crit200", log_x = True, log_y = True)
    #produce_simulations_graph(stellar_mass, y_axis_label = "M* ($M_{sun}$)", graph_title_partial = "Central Subhalo Stellar Mass", log_x = True, log_y = True)
    #produce_simulations_graph(black_hole_mass, y_axis_label = "$M_{BH}$ ($M_{sun}$)", graph_title_partial = "Central Subhalo Black Hole Mass", log_x = True, log_y = True)
    #produce_simulations_graph(dark_matter_mass, y_axis_label = "$M_{DM}$ ($M_{sun}$)", graph_title_partial = "Central Subhalo Dark Matter Mass", log_x = True, log_y = True)
    #produce_simulations_graph(baryon_mass, y_axis_label = "$M_{Baryonic}$ ($M_{sun}$)", graph_title_partial = "Central Subhalo Baryonic Mass", log_x = True, log_y = True)
    #t = -7#8
    #print(assembily_history[Simulations.Early][constants.tags[t]])
    #a = gas_mass(assembily_history[Simulations.Early][constants.tags[t]]["halo"], assembily_history[Simulations.Early][constants.tags[t]]["subhalo"], constants.tags[t], Simulations.Early)
    produce_single_simulation_graphs([dark_matter_mass, gas_mass, stellar_mass, black_hole_mass], ["Dark Matter", "Gas", "Stars", "Black Holes"], y_axis_label = "Mass ($M_{sun}$)", graph_title_partial = "Central Subhalo Mass Brakedown (r = 30 kPc)", log_x = True, log_y = True)

    #produce_single_simulation_graphs([total_mass_M200, dark_matter_mass, baryon_mass, lambda *args, **kwargs: dark_matter_mass(*args, **kwargs) + baryon_mass(*args, **kwargs)], ["Halo $M_{200}$", "Dark Matter (r = 30 kPc)", "Baryonic Matter (r = 30 kPc)", "All Matter (r = 30 kPc)"], y_axis_label = "Mass ($M_{sun}$)", graph_title_partial = "Central Subhalo Mass", log_x = True, log_y = True)
