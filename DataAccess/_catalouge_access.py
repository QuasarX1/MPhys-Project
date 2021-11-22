import numpy as np
import h5py as h5
import os

from ._unit_conversions import UnitSystem
from ._simulation_combinations import Simulations, SimulationModels

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
    for file_number in file_numbers:
        with h5.File(subfind_file_template.format(file_number), "r") as datafile:
            data = datafile[f"/{table}/{field_name}"]

            if file_number == 0:
                h_scale_exponent = data.attrs["h-scale-exponent"]
                a_scale_exponent = data.attrs["aexp-scale-exponent"]
                cgs_conversion_factor = data.attrs["CGSConversionFactor"]
                h = datafile["Header"].attrs["HubbleParam"]
                a = datafile["Header"].attrs["ExpansionFactor"]

                if verbose:
                    print("Loading", field_name)
                    print("h exponent =", h_scale_exponent)
                    print("a exponent =", a_scale_exponent)
                    print("cgs conversion factor =", cgs_conversion_factor)
                    print("h =", h)
                    print("a =", a)
                    
                dt = data.dtype
                data_arr = np.array(data, dtype = dt)
            
            else:
                data_arr = np.append(data_arr,np.array(data, dtype = dt), axis = 0)

    if verbose: print("Run out of files after loading", file_number)

    # Convert values - no corrections needed for integer type data
    if not np.issubdtype(dt, np.integer):
        # cgs numbers can be huge and overflow np.float32
        # Recast the data to float64 to be safe
        if unit_system == UnitSystem.cgs:
            data_arr = np.array(data_arr, dtype = np.float64)
        data_arr = UnitSystem.convert_data(data_arr, UnitSystem.h_less_comoving_GADGET, unit_system, h, h_scale_exponent, a, a_scale_exponent, cgs_conversion_factor)

        #if unit_system == UnitSystem.physical:
        #    data_arr *= np.power(h, h_scale_exponent) * np.power(a, a_scale_exponent)
        #
        #if unit_system == UnitSystem.cgs:
        #
        #    # cgs numbers can be huge and overflow np.float32
        #    # Recast the data to float64 to be safe
        #    data_arr = np.array(data_arr, dtype = np.float64)
        #
        #    data_arr *= cgs_conversion_factor

    return data_arr
