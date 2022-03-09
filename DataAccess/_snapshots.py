import numpy as np
import pyread_eagle
import h5py as h5
import os

from DataAccess._unit_conversions import UnitSystem
from DataAccess._simulation_combinations import Simulations, SimulationModels
from DataAccess._particle_type import ParticleType

def combine_limits(lower, upper):
    bound_length = len(lower)
    if bound_length != len(upper):
        raise ValueError("Arrays of bounds were not of the same length.")
    return np.append(np.array(lower).reshape((bound_length, 1)), np.array(upper).reshape((bound_length, 1)), axis = 1).reshape(2 * bound_length)

class ParticleReadConversion_EagleSnapshot(pyread_eagle.EagleSnapshot):
    """
    Custom Eagle Smapshot class that implements extra methods for reading data and doing unit conversions.

    Constructor Paramiters:
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
                      particle_type:         ParticleType,
                      field_name:            str,
                      lower_limits:          np.ndarray   = None,
                      upper_limits:          np.ndarray   = None,
                      unit_system:           UnitSystem   = UnitSystem.h_less_comoving_GADGET,
                      limits_unit_system:    UnitSystem   = UnitSystem.h_less_comoving_GADGET,
                      use_vanilla_selection: bool         = False,
                      return_coordinates:    bool         = False) -> np.ndarray:
        """
        Loads particle data in the set region for a specified snapshot.

        See https://j-davies-ari.github.io/eagle-guide/working_with_particles for origanal code and tutorial.

        Paramiters:
            ParticleType particle_type         -> The type of particle for which to load data
                     str field_name            -> The name of the data field
              np.ndarray lower_limits          -> Lower boundries of the region to select particles within (an array with 3 values - [x, y, z]) (default: lower limits of the box)
              np.ndarray upper_limits          -> Upper boundries of the region to select particles within (an array with 3 values - [x, y, z]) (default: upper limits of the box)
              UnitSystem unit_system           -> The unit system to convert the results into (default: UnitSystem.h_less_comoving_GADGET)
              UnitSystem limits_unit_system    -> The unit system that represents the values of the specified limits (default: UnitSystem.h_less_comoving_GADGET)
                    bool use_vanilla_selection -> Ignore specified limits and fall back on pyread_eagle.EagleSnapshot selection (default: False)
                    bool return_coordinates    -> Return the coordinates associated with each particle (default: False)

        Returns:
            np.ndarray -> Numpy array containing the data from the specified field
            np.ndarray -> Optionally provides an array of coordinates
        """

        if not use_vanilla_selection:
            if lower_limits is None:
                # If no limits are specified, use the box boundry
                lower_limits = np.full((3, ), 0.0)
            elif limits_unit_system != UnitSystem.h_less_comoving_GADGET:
                # If limits are specified, ensure they are in GADGET units
                lower_limits = self.convert_distance_values(lower_limits, limits_unit_system, UnitSystem.h_less_comoving_GADGET)

            if upper_limits is None:
                # If no limits are specified, use the box boundry
                upper_limits = np.full((3, ), self.boxsize)
            elif limits_unit_system != UnitSystem.h_less_comoving_GADGET:
                # If limits are specified, ensure they are in GADGET units
                upper_limits = self.convert_distance_values(upper_limits, limits_unit_system, UnitSystem.h_less_comoving_GADGET)

            limits = combine_limits(lower_limits, upper_limits)
            self.select_region(*limits)

        h_scale_exponent = self.h_scale_exponent_values[particle_type][field_name]
        a_scale_exponent = self.a_scale_exponent_values[particle_type][field_name]
        cgs_conversion_factor = self.cgs_conversion_factor_values[particle_type][field_name]

        # Load the data
        data_arr = self.read_dataset(particle_type.value, field_name)
        particle_coordiantes = self.read_dataset(particle_type.value, "Coordinates")
        # Convert to numpy array with the correct data type
        dt = data_arr.dtype
        data_arr = np.array(data_arr, dtype = dt)
        if not use_vanilla_selection:
            # Remove excess particles that were loaded
            filter = (particle_coordiantes[:, 0] > lower_limits[0]) & (particle_coordiantes[:, 0] < upper_limits[0]) & (particle_coordiantes[:, 1] > lower_limits[1]) & (particle_coordiantes[:, 1] < upper_limits[1]) & (particle_coordiantes[:, 2] > lower_limits[2]) & (particle_coordiantes[:, 2] < upper_limits[2])
            data_arr = data_arr[filter]
            particle_coordiantes = particle_coordiantes[filter]


        # Convert values - no corrections needed for integer type data
        if not np.issubdtype(dt, np.integer):
            # cgs numbers can be huge and overflow np.float32
            # Recast the data to float64 to be safe
            if unit_system == UnitSystem.cgs:
                data_arr = np.array(data_arr, dtype = np.float64)
            data_arr = UnitSystem.convert_data(data_arr, UnitSystem.h_less_comoving_GADGET, unit_system, self.hubble_paramiter, h_scale_exponent, self.expansion_factor, a_scale_exponent, cgs_conversion_factor)

        # Coordinates will require conversion
        if return_coordinates:
            # cgs numbers can be huge and overflow np.float32
            # Recast the data to float64 to be safe
            if unit_system == UnitSystem.cgs:
                particle_coordiantes = np.array(particle_coordiantes, dtype = np.float64)
            particle_coordiantes = UnitSystem.convert_data(particle_coordiantes, UnitSystem.h_less_comoving_GADGET, unit_system, self.hubble_paramiter, self.h_scale_exponent_values[particle_type]["Coordinates"], self.expansion_factor, self.a_scale_exponent_values[particle_type]["Coordinates"], self.cgs_conversion_factor_values[particle_type]["Coordinates"])


        if not return_coordinates:
            return data_arr
        else:
            return data_arr, particle_coordiantes

    def centre_particles(self, coordinates, centre):
        """
        Changes the raw coordinates provided so they are relitive to a central point and the periodic box accounted for.

        Paramiters:
            np.ndarray coordinates -> Array of coordinate arrays (an array with 3 values - [x, y, z]) to be re-expressed
            np.ndarray centre      -> Point to express coordinates relitive to (an array with 3 values - [x, y, z])

        Returns:
            np.ndarray -> New array of coordinate arrays (an array with 3 values - [x, y, z]) expressed relitive to the centre point
        """
        return ((coordinates - centre + (self.boxsize / 2)) % self.boxsize) - (self.boxsize / 2)

    def particle_read_sphere(self,
                             particle_type:         ParticleType,
                             field_name:            str,
                             centre:                np.ndarray,
                             radius:                float,
                             unit_system:           UnitSystem   = UnitSystem.h_less_comoving_GADGET,
                             limits_unit_system:    UnitSystem   = UnitSystem.h_less_comoving_GADGET,
                             return_coordinates:    bool         = False) -> np.ndarray:
        """
        Loads particle data in the set region for a specified snapshot.

        See https://j-davies-ari.github.io/eagle-guide/working_with_particles for origanal code and tutorial.

        Paramiters:
            ParticleType particle_type         -> The type of particle for which to load data
                     str field_name            -> The name of the data field
              np.ndarray centre                -> Centre of the region to select particles within (an array with 3 values - [x, y, z])
                   float radius                -> Radius of the region to select particles within
              UnitSystem unit_system           -> The unit system to convert the results into (default: UnitSystem.h_less_comoving_GADGET)
              UnitSystem limits_unit_system    -> The unit system that represents the values of the specified limits (default: UnitSystem.h_less_comoving_GADGET)
                    bool use_vanilla_selection -> Ignore specified limits and fall back on pyread_eagle.EagleSnapshot selection (default: False)
                    bool return_coordinates    -> Return the coordinates associated with each particle (default: False)

        Returns:
            np.ndarray -> Numpy array containing the data from the specified field
            np.ndarray -> Optionally provides an array of coordinates
        """

        if limits_unit_system != UnitSystem.h_less_comoving_GADGET:
            # Ensure limits are in GADGET units
            centre = self.convert_distance_values(centre, limits_unit_system, UnitSystem.h_less_comoving_GADGET)
            radius = self.convert_distance_values(radius, limits_unit_system, UnitSystem.h_less_comoving_GADGET)

        limits = combine_limits(centre - radius, centre + radius)
        self.select_region(*limits)

        h_scale_exponent = self.h_scale_exponent_values[particle_type][field_name]
        a_scale_exponent = self.a_scale_exponent_values[particle_type][field_name]
        cgs_conversion_factor = self.cgs_conversion_factor_values[particle_type][field_name]

        # Load the data
        data_arr = self.read_dataset(particle_type.value, field_name)
        coordinates = self.read_dataset(particle_type.value, "Coordinates")

        # Convert to numpy array with the correct data type
        dt = data_arr.dtype
        data_arr = np.array(data_arr, dtype = dt)
        # Centre the region coordinates
        coordinates = self.centre_particles(coordinates, centre)
        # Remove excess particles that were loaded
        data_arr = data_arr[(coordinates**2).sum(axis = 1) < radius**2]
        coordinates = coordinates[(coordinates**2).sum(axis = 1) < radius**2]


        # Convert values - no corrections needed for integer type data
        if not np.issubdtype(dt, np.integer):
            # cgs numbers can be huge and overflow np.float32
            # Recast the data to float64 to be safe
            if unit_system == UnitSystem.cgs:
                data_arr = np.array(data_arr, dtype = np.float64)
            data_arr = UnitSystem.convert_data(data_arr, UnitSystem.h_less_comoving_GADGET, unit_system, self.hubble_paramiter, h_scale_exponent, self.expansion_factor, a_scale_exponent, cgs_conversion_factor)

        # Coordinates will require conversion
        if return_coordinates:
            # cgs numbers can be huge and overflow np.float32
            # Recast the data to float64 to be safe
            if unit_system == UnitSystem.cgs:
                coordinates = np.array(coordinates, dtype = np.float64)
            coordinates = UnitSystem.convert_data(coordinates, UnitSystem.h_less_comoving_GADGET, unit_system, self.hubble_paramiter, self.h_scale_exponent_values[particle_type]["Coordinates"], self.expansion_factor, self.a_scale_exponent_values[particle_type]["Coordinates"], self.cgs_conversion_factor_values[particle_type]["Coordinates"])


        if not return_coordinates:
            return data_arr
        else:
            return data_arr, coordinates

    def convert_particle_values(self, data, from_unit, to_unit, particle_type, field_name):
        return UnitSystem.convert_data(data, from_unit, to_unit, self.hubble_paramiter, self.h_scale_exponent_values[particle_type][field_name], self.expansion_factor, self.a_scale_exponent_values[particle_type][field_name], self.cgs_conversion_factor_values[particle_type][field_name])

    def convert_distance_values(self, data, from_unit, to_unit):
        #TODO: check following assumption
        # These shoould have the same conversion factors as distance
        particle_type = ParticleType.gas
        field_name = "Coordinates"

        return UnitSystem.convert_data(data, from_unit, to_unit, self.hubble_paramiter, self.h_scale_exponent_values[particle_type][field_name], self.expansion_factor, self.a_scale_exponent_values[particle_type][field_name], self.cgs_conversion_factor_values[particle_type][field_name])
