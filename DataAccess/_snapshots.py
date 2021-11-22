import numpy as np
import pyread_eagle
import h5py as h5
import os

from ._unit_conversions import UnitSystem
from ._simulation_combinations import Simulations, SimulationModels
from ._particle_type import ParticleType

class ParticleReadConversion_EagleSnapshot(pyread_eagle.EagleSnapshot):
    """
    Custom Eagle Smapshot class that implements a method
                     str tag        -> Snapshot identifier
             Simulations simulation -> The target simulation specified with an enum
        SimulationModels model      -> The target simulation model specified with an enum
                     str data_root  -> The filepath to the root of the GM simulations directory
    """

    def __init__(self,
                 tag:        str,
                 simulation: Simulations,
                 model:      SimulationModels,
                 data_root:  str):
        self.__snapshot_folder = os.path.join(data_root, SimulationModels.to_string(model), Simulations.to_string(simulation), f"snapshot_{tag}")
        self.__snapfile_template = os.path.join(self.__snapshot_folder, f"snap_{tag}.{{}}.hdf5")
        self.__snapfile = self.__snapfile_template.format(0)
        super().__init__(self.__snapfile)

    def particle_read(self,
                      particle_type: ParticleType,
                      field_name:    str,
                      unit_system:   UnitSystem = UnitSystem.h_less_comoving_GADGET,
                      verbose:       bool       = False) -> np.ndarray:
        """
        Loads particle data in the set region for a specified snapshot.

        See https://j-davies-ari.github.io/eagle-guide/working_with_particles for origanal code and tutorial.

        Paramiters:
            ParticleType particle_type -> The type of particle for which to load data
                     str field_name    -> The name of the data field
              UnitSystem unit_system   -> The unit system to convert the results into (default: UnitSystem.h_less_comoving_GADGET)
                    bool verbose       -> Print out progress infomation (default: False)

        Returns:
            np.ndarray -> Numpy array containing the data from the specified field
        """

        particle_type_number = particle_type.value
        
        file_numbers = [int(file_name[18:-5]) for file_name in os.listdir(self.__snapshot_folder) if file_name[:5] == "snap_"]
        file_numbers.sort()
        dt = None
        data_arr = None
        h_scale_exponent = None
        a_scale_exponent = None
        cgs_conversion_factor = None
        h = None
        a = None
        for file_number in file_numbers:
            with h5.File(self.__snapfile_template.format(file_number), "r") as datafile:
                try:
                    # Read file header values to allow conversions
                    h = datafile["Header"].attrs["HubbleParam"]
                    a = datafile["Header"].attrs["ExpansionFactor"]
                    h_scale_exponent = datafile[f"/PartType{particle_type_number}/{field_name}"].attrs["h-scale-exponent"]
                    a_scale_exponent = datafile[f"/PartType{particle_type_number}/{field_name}"].attrs["aexp-scale-exponent"]
                    cgs_conversion_factor = datafile[f"/PartType{particle_type_number}/{field_name}"].attrs["CGSConversionFactor"]

                    if verbose:
                        print("Loading", field_name)
                        print("h exponent =", h_scale_exponent)
                        print("a exponent =", a_scale_exponent)
                        print("cgs conversion factor =", cgs_conversion_factor)
                        print("h =", h)
                        print("a =", a)

                    break

                except:
                    # If no particles of the type specified are present, try the next file header for the values
                    if verbose: print("No particles of type", ptype, "in snapfile", snapfile_ind)

        # Load the data
        data_arr = self.read_dataset(particle_type.value, field_name)
        # Convert to numpy array with the correct data type
        dt = data_arr.dtype
        data_arr = np.array(data_arr, dtype = dt)


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
