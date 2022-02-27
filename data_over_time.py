from enum import Enum
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




class X_Axis_Value(Enum):
    redshift = 0
    time = 1
    expansion_factor = 2

def produce_simulations_graph(func, y_axis_label, graph_title_partial, x_axis: X_Axis_Value = X_Axis_Value.redshift, log_x = True, log_y = False, invert_x = True, invert_y = False, xlim_overide = (None, None), ylim_overide = (None, None)):
    quantity_values = []
    redshift_values = []
    time_values = []
    expansion_factor_values = []
    for i, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
        quantity_values.append(np.empty(shape = (len(constants.tags), )))
        redshift_values.append(np.array([float("{}.{}".format(tag[5:8], tag[9:12])) for tag in constants.tags]))
        time_values.append(np.array(list(constants.times.values())))
        expansion_factor_values.append(np.array(list(constants.expansion_factors.values())))
        for j, tag in enumerate(constants.tags):
            try:
                halo = assembily_history[simulation][tag]["halo"]
                subhalo = assembily_history[simulation][tag]["subhalo"]
                if halo is None:
                    raise LookupError("No history avalible.")

                quantity_values[i][j] = func(halo, subhalo, tag, simulation)

            except LookupError:
                quantity_values[i][j] = None
            
        array_filter = True ^ np.isnan(quantity_values[i])
        quantity_values[i] = quantity_values[i][array_filter]
        redshift_values[i] = redshift_values[i][array_filter]
        time_values[i] = time_values[i][array_filter]
        expansion_factor_values[i] = expansion_factor_values[i][array_filter]

        
    # Set value for x-axis and label
    x_values = None
    if x_axis == X_Axis_Value.redshift:
        x_values = redshift_values
        plt.xlabel("Redshift")
    elif x_axis == X_Axis_Value.time:
        x_values = time_values
        plt.xlabel("Time (Gyrs)")
    elif x_axis == X_Axis_Value.expansion_factor:
        x_values = expansion_factor_values
        plt.xlabel("Expansion Factor")
    else:
        raise ValueError("Value of \"x_axis\" was not from the \"X_Axis_Value\" enum.")

    plt.plot(x_values[0], quantity_values[0], label = "Early")
    plt.plot(x_values[1], quantity_values[1], label = "Organic")
    plt.plot(x_values[2], quantity_values[2], label = "Late")

    # Set axis log status
    if log_x and log_y:
        plt.loglog()
    elif log_x:
        plt.semilogx()
    elif log_y:
        plt.semilogy()

    # Set axis limits
    if xlim_overide[0] is not None:
        plt.xlim(xlim_overide[0], plt.xlim()[1])
    if xlim_overide[1] is not None:
        plt.xlim(plt.xlim()[0], xlim_overide[1])

    if ylim_overide[0] is not None:
        plt.ylim(ylim_overide[0], plt.ylim()[1])
    if ylim_overide[1] is not None:
        plt.ylim(plt.ylim()[0], ylim_overide[1])

    # Invert axes
    if invert_x:
        plt.xlim(plt.xlim()[1], plt.xlim()[0])
    if invert_y:
        plt.ylim(plt.ylim()[1], plt.ylim()[0])
        
    plt.ylabel(y_axis_label)
    plt.title(f"{graph_title_partial} over Time")
    plt.legend()
    plt.show()

def produce_single_simulation_graphs(funcs, quantity_labels, y_axis_label, graph_title_partial, x_axis: X_Axis_Value = X_Axis_Value.redshift, log_x = True, log_y = False, invert_x = True, invert_y = False, xlim_overide = (None, None), ylim_overide = (None, None)):
    fig, axes = plt.subplots(1, 3)

    for simulation_number, simulation in enumerate((Simulations.Early, Simulations.Organic, Simulations.Late)):
        quantity_values = []
        redshift_values = []
        time_values = []
        expansion_factor_values = []
        for i in range(len(funcs)):
            quantity_values.append(np.empty(shape = (len(constants.tags), )))
            redshift_values.append(np.array([float("{}.{}".format(tag[5:8], tag[9:12])) for tag in constants.tags]))
            time_values.append(np.array(list(constants.times.values())))
            expansion_factor_values.append(np.array(list(constants.expansion_factors.values())))
        
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
            redshift_values[i] = redshift_values[i][array_filter]
            time_values[i] = time_values[i][array_filter]
            expansion_factor_values[i] = expansion_factor_values[i][array_filter]

            # Set value for x-axis and label
            x_values = None
            if x_axis == X_Axis_Value.redshift:
                x_values = redshift_values
                axes[simulation_number].set_xlabel("Redshift")
            elif x_axis == X_Axis_Value.time:
                x_values = time_values
                axes[simulation_number].set_xlabel("Time (Gyrs)")
            elif x_axis == X_Axis_Value.expansion_factor:
                x_values = expansion_factor_values
                axes[simulation_number].set_xlabel("Expansion Factor")
            else:
                raise ValueError("Value of \"x_axis\" was not from the \"X_Axis_Value\" enum.")
        
            axes[simulation_number].plot(x_values[i], quantity_values[i], label = quantity_labels[i])

        # Set axis log status
        if log_x and log_y:
            axes[simulation_number].loglog()
        elif log_x:
            axes[simulation_number].semilogx()
        elif log_y:
            axes[simulation_number].semilogy()

        # Set axis limits
        if xlim_overide[0] is not None:
            axes[simulation_number].set_xlim(xlim_overide[0], axes[simulation_number].get_xlim()[1])
        if xlim_overide[1] is not None:
            axes[simulation_number].set_xlim(axes[simulation_number].get_xlim()[0], xlim_overide[1])

        if ylim_overide[0] is not None:
            axes[simulation_number].set_ylim(ylim_overide[0], axes[simulation_number].get_ylim()[1])
        if ylim_overide[1] is not None:
            axes[simulation_number].set_ylim(axes[simulation_number].get_ylim()[0], ylim_overide[1])

        # Invert axes
        if invert_x:
            axes[simulation_number].set_xlim(axes[simulation_number].get_xlim()[1], axes[simulation_number].get_xlim()[0])
        if invert_y:
            axes[simulation_number].set_ylim(axes[simulation_number].get_ylim()[1], axes[simulation_number].get_ylim()[0])

        axes[simulation_number].set_ylabel(y_axis_label)
        simulation_name = Simulations.to_string(simulation)
        #axes[simulation_number].title(f"{simulation_name} {graph_title_partial} for each Snapshot")
        axes[simulation_number].set_title(f"{simulation_name}")
        axes[simulation_number].legend()
    plt.suptitle(f"{graph_title_partial} over Time")
    plt.show()



if __name__ == "__main__":
    # M_200 (r = R_200)
    produce_simulations_graph(total_mass_M200, y_axis_label = "$M_{200}$ ($M_{sun}$)", graph_title_partial = "Group_M_Crit200",
                              x_axis = X_Axis_Value.time, log_x = False, log_y = True, invert_x = False, ylim_overide = (10**10, None))
    
    # M* (r = 30KPc)
    produce_simulations_graph(stellar_mass, y_axis_label = "$M_*$ ($M_{sun}$)", graph_title_partial = "Central Subhalo Stellar Mass",
                              x_axis = X_Axis_Value.time, log_x = False, log_y = True, invert_x = False, ylim_overide = (10**6, None))
    
    # M_BH (r = 30KPc)
    produce_simulations_graph(black_hole_mass, y_axis_label = "$M_{BH}$ ($M_{sun}$)", graph_title_partial = "Central Subhalo Black Hole Mass",
                              x_axis = X_Axis_Value.time, log_x = False, log_y = True, invert_x = False, ylim_overide = (10**5, None))
    
    # M_DM (r = 30KPc)
    produce_simulations_graph(dark_matter_mass, y_axis_label = "$M_{DM}$ ($M_{sun}$)", graph_title_partial = "Central Subhalo Dark Matter Mass",
                              x_axis = X_Axis_Value.time, log_x = False, log_y = True, invert_x = False)
    
    # Barionic Mass (r = 30KPc)
    produce_simulations_graph(baryon_mass, y_axis_label = "$M_{Baryonic}$ ($M_{sun}$)", graph_title_partial = "Central Subhalo Baryonic Mass",
                              x_axis = X_Axis_Value.time, log_x = False, log_y = True, invert_x = False)
    
    # Mass Brakedown (r = 30KPc)
    produce_single_simulation_graphs([dark_matter_mass, gas_mass, stellar_mass, black_hole_mass], ["Dark Matter", "Gas", "Stars", "Black Holes"], y_axis_label = "Mass ($M_{sun}$)", graph_title_partial = "Central Subhalo Mass Brakedown (r = 30 kPc)",
                                     x_axis = X_Axis_Value.time, log_x = False, log_y = True, invert_x = False)
