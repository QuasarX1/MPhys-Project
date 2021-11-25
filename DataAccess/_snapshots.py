import numpy as np
import pyread_eagle
import h5py as h5
import os

from DataAccess._unit_conversions import UnitSystem
from DataAccess._simulation_combinations import Simulations, SimulationModels
from DataAccess._particle_type import ParticleType

class ParticleReadConversion_EagleSnapshot(pyread_eagle.EagleSnapshot):
    """
    Custom Eagle Smapshot class that implements a method
                     str tag        -> Snapshot identifier
             Simulations simulation -> The target simulation specified with an enum
        SimulationModels model      -> The target simulation model specified with an enum
                     str data_root  -> The filepath to the root of the GM simulations directory
                    bool verbose    -> Print out progress infomation (default: False)
    """

    def __init__(self,
                 tag:        str,
                 simulation: Simulations,
                 model:      SimulationModels,
                 data_root:  str,
                 verbose:    bool             = False):
        self.__snapshot_folder = os.path.join(data_root, SimulationModels.to_string(model), Simulations.to_string(simulation), f"snapshot_{tag}")
        self.__snapfile_template = os.path.join(self.__snapshot_folder, f"snap_{tag}.{{}}.hdf5")
        self.__snapfile = self.__snapfile_template.format(0)
        with h5.File(self.__snapfile, 'r') as data_file:
            self.header = dict(data_file['/Header'].attrs)
        super().__init__(self.__snapfile, verbose)

        particle_types = (ParticleType.gas, ParticleType.dark_matter, ParticleType.star, ParticleType.black_hole)

        self.h_scale_exponent_values = {particle_type: {} for particle_type in particle_types}
        self.a_scale_exponent_values = {particle_type: {} for particle_type in particle_types}
        self.cgs_conversion_factor_values = {particle_type: {} for particle_type in particle_types}
        
        file_numbers = [int(file_name[18:-5]) for file_name in os.listdir(self.__snapshot_folder) if file_name[:5] == "snap_"]
        file_numbers.sort()
        for particle_type in particle_types:
            particle_type_number = particle_type.value
            for file_number in file_numbers:
                with h5.File(self.__snapfile_template.format(file_number), "r") as datafile:
                    try:
                        # Read file header values to allow conversions
                        for field_name in dict(datafile[f"/PartType{particle_type_number}"]):
                            # Ignore this field if it dosen't have conversion exponents or factors
                            if "h-scale-exponent" not in datafile[f"/PartType{particle_type_number}/{field_name}"].attrs:
                                continue

                            self.h_scale_exponent_values[particle_type][field_name] = datafile[f"/PartType{particle_type_number}/{field_name}"].attrs["h-scale-exponent"]
                            self.a_scale_exponent_values[particle_type][field_name] = datafile[f"/PartType{particle_type_number}/{field_name}"].attrs["aexp-scale-exponent"]
                            self.cgs_conversion_factor_values[particle_type][field_name] = datafile[f"/PartType{particle_type_number}/{field_name}"].attrs["CGSConversionFactor"]

                        break

                    except:# is this a KeyError? - should occour ONLY when the "field" is not present in the "PartType{n}" group as valid fields without needed attributes are handled?
                        # If no particles of the type specified are present, try the next file header for the values
                        if verbose: print("No particles of type", particle_type_number, "in snapfile", file_number)

    @property
    def hubble_paramiter(self):
        return self.header['HubbleParam']

    @property
    def expansion_factor(self):
        return self.header['ExpansionFactor']

    def particle_read(self,
                      particle_type: ParticleType,
                      field_name:    str,
                      unit_system:   UnitSystem = UnitSystem.h_less_comoving_GADGET) -> np.ndarray:
        """
        Loads particle data in the set region for a specified snapshot.

        See https://j-davies-ari.github.io/eagle-guide/working_with_particles for origanal code and tutorial.

        Paramiters:
            ParticleType particle_type -> The type of particle for which to load data
                     str field_name    -> The name of the data field
              UnitSystem unit_system   -> The unit system to convert the results into (default: UnitSystem.h_less_comoving_GADGET)

        Returns:
            np.ndarray -> Numpy array containing the data from the specified field
        """

        h_scale_exponent = self.h_scale_exponent_values[particle_type][field_name]
        a_scale_exponent = self.a_scale_exponent_values[particle_type][field_name]
        cgs_conversion_factor = self.cgs_conversion_factor_values[particle_type][field_name]

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
            data_arr = UnitSystem.convert_data(data_arr, UnitSystem.h_less_comoving_GADGET, unit_system, self.hubble_paramiter, h_scale_exponent, self.expansion_factor, a_scale_exponent, cgs_conversion_factor)

        return data_arr

    def convert_particle_values(self, data, from_unit, to_unit, particle_type, field_name):
        return UnitSystem.convert_data(data, from_unit, to_unit, self.hubble_paramiter, self.h_scale_exponent_values[particle_type][field_name], self.expansion_factor, self.a_scale_exponent_values[particle_type][field_name], self.cgs_conversion_factor_values[particle_type][field_name])

    def convert_distance_values(self, data, from_unit, to_unit):
        #TODO: check following assumption
        # These shoould have the same conversion factors as distance
        particle_type = ParticleType.gas
        field_name = "Coordinates"

        return UnitSystem.convert_data(data, from_unit, to_unit, self.hubble_paramiter, self.h_scale_exponent_values[particle_type][field_name], self.expansion_factor, self.a_scale_exponent_values[particle_type][field_name], self.cgs_conversion_factor_values[particle_type][field_name])
