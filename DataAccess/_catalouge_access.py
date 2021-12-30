import numpy as np
import h5py as h5
import os

from DataAccess._unit_conversions import UnitSystem
from DataAccess._particle_type import ParticleType
from DataAccess._simulation_combinations import Simulations, SimulationModels

def load_catalouge_field(field_name:  str,
                         table:       str,
                         tag:         str,
                         simulation:  Simulations,
                         model:       SimulationModels,
                         data_root:   str,
                         unit_system: UnitSystem       = UnitSystem.h_less_comoving_GADGET,
                         verbose:     bool             = False) -> np.ndarray:
    """
    Loads the data from a specified field for a specified snapshot.

    NOTE: only "FOF" and "Subhalo" tables are supported.

    See https://j-davies-ari.github.io/eagle-guide/working_with_catalogues for origanal code and tutorial.

    Paramiters:
                     str field_name  -> The name of the data field
                     str table       -> The name of the table
                     str tag         -> Snapshot identifier
             Simulations simulation  -> The target simulation specified with an enum
        SimulationModels model       -> The target simulation model specified with an enum
                     str data_root   -> The filepath to the root of the GM simulations directory
              UnitSystem unit_system -> The unit system to convert the results into (default: UnitSystem.h_less_comoving_GADGET)
                    bool verbose     -> Print out progress infomation (default: False)

    Returns:
        np.ndarray -> Numpy array containing the data from the specified field
    """

    assert table in ["FOF","Subhalo"], "table must be either FOF or Subhalo"

    

    group_folder = os.path.join(data_root, SimulationModels.to_string(model), Simulations.to_string(simulation), f"groups_{tag}")
    subfind_file_template = os.path.join(group_folder, f"eagle_subfind_tab_{tag}.{{}}.hdf5")

    file_numbers = [int(file_name[31:-5]) for file_name in os.listdir(group_folder) if file_name[:17] == "eagle_subfind_tab"]
    file_numbers.sort()
    dt = None
    data_arr = None
    h_scale_exponent = None
    a_scale_exponent = None
    cgs_conversion_factor = None
    h = None
    a = None
    read_constants = True
    for file_number in file_numbers:
        with h5.File(subfind_file_template.format(file_number), "r") as datafile:
            try:
                data = datafile[f"/{table}/{field_name}"]

                #if file_number == 0:
                if read_constants:
                    try:
                        h_scale_exponent = data.attrs["h-scale-exponent"]
                        a_scale_exponent = data.attrs["aexp-scale-exponent"]
                        cgs_conversion_factor = data.attrs["CGSConversionFactor"]
                        h = datafile["Header"].attrs["HubbleParam"]
                        a = datafile["Header"].attrs["ExpansionFactor"]
                        read_constants = False

                        if verbose:
                            print("Loading", field_name)
                            print("h exponent =", h_scale_exponent)
                            print("a exponent =", a_scale_exponent)
                            print("cgs conversion factor =", cgs_conversion_factor)
                            print("h =", h)
                            print("a =", a)

                    except:
                        read_constants = False

                        if verbose:
                            print("Error reading data scaling constants! Will attempt again with next file.")
                    
                    dt = data.dtype
                    data_arr = np.array(data, dtype = dt)
            
                else:
                    data_arr = np.append(data_arr,np.array(data, dtype = dt), axis = 0)

            except KeyError:
                if verbose:
                    print(f"No data avalible fron file number {file_number}.")

    if verbose: print("Run out of files after loading", file_number)

    if read_constants:# No files contained the nessessary field or the constants failed to load
        if verbose:
            print("Unable to retrive scaling constants - does this tag contain the specified field?")
        raise LookupError("Scaling constants were not loaded from tag files. Possible that the specified field was not present.")

    # Convert values - no corrections needed for integer type data
    if not np.issubdtype(dt, np.integer):
        # cgs numbers can be huge and overflow np.float32
        # Recast the data to float64 to be safe
        if unit_system == UnitSystem.cgs:
            data_arr = np.array(data_arr, dtype = np.float64)
        data_arr = UnitSystem.convert_data(data_arr, UnitSystem.h_less_comoving_GADGET, unit_system, h, h_scale_exponent, a, a_scale_exponent, cgs_conversion_factor)

    return data_arr

def particles_by_group(tag:           str,
                       particle_type: ParticleType,
                       simulation:    Simulations,
                       model:         SimulationModels,
                       data_root:     str,
                       verbose:       bool             = False):
    """
    Load the particle IDs, subgroup number and group number and the binding energy for particles in a tag.

    Paramiters:
                     str tag           -> Snapshot identifier
            ParticleType particle_type -> The type of particle to load data for, specified with an enum
             Simulations simulation    -> The target simulation specified with an enum
        SimulationModels model         -> The target simulation model specified with an enum
                     str data_root     -> The filepath to the root of the GM simulations directory
                    bool verbose       -> Print out progress infomation (default: False)

    Returns:
        np.ndarray -> Numpy array containing the particle ID, the associated group number, the associated subgroup number and the binding energy
    """

    subfind_particles_folder = os.path.join(data_root, SimulationModels.to_string(model), Simulations.to_string(simulation), f"particledata_{tag}")
    subfind_file_template = os.path.join(subfind_particles_folder, f"eagle_subfind_particles_{tag}.{{}}.hdf5")

    file_numbers = [int(file_name[37:-5]) for file_name in os.listdir(subfind_particles_folder) if file_name[:23] == "eagle_subfind_particles"]
    file_numbers.sort()
    results = [[], [], [], []]
    for file_number in file_numbers:
        with h5.File(subfind_file_template.format(file_number), "r") as datafile:
            particle_type_number = particle_type.value
            try:
                results[0].append(datafile[f"/PartType{particle_type_number}/ParticleIDs"])# Particles
                results[1].append(datafile[f"/PartType{particle_type_number}/GroupNumber"])# Associated Groups
                results[2].append(datafile[f"/PartType{particle_type_number}/SubGroupNumber"])# Associated Subgroups
                results[3].append(datafile[f"/PartType{particle_type_number}/ParticleBindingEnergy"])# Associated Binding Energy
            except KeyError:
                if verbose:
                    print(f"No data avalible fron file number {file_number}.")

            if verbose: print("Run out of files after loading", file_number)

    return np.array(results)