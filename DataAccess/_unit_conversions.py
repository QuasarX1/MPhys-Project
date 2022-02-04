from enum import Enum
import numpy as np

from DataAccess._constants import SimulationConstants

class UnitSystem(Enum):
    """
    Systems of units that are accepted paramiters for functions in the "DataAccess" package.
    """
    h_less_comoving_GADGET = 0
    physical = 1
    cgs = 2
    physical_comoving = 3
    cgs_comoving = 4

    @staticmethod
    def convert_data(data: np.ndarray, from_unit, to_unit, h = None, h_scale_exponent = None, a = None, a_scale_exponent = None, cgs_conversion_factor = None) -> np.ndarray:
        """
        Converts data between two different supported unit systems.

        Paramiters:
            (np.ndarray or float) data                  -> The value(s) to be converted
                       UnitSystem from_unit             -> The system of units by which the data is currently described
                       UnitSystem to_unit               -> The target system of units
                            float h                     -> Dimensionless Hubble Paramiter
                            float h_scale_exponent      -> Exponent to apply to the value of h
                            float a                     -> Expansion Factor
                            float a_scale_exponent      -> Exponent to apply to the value of a
                            float cgs_conversion_factor -> Conversion factor between physical and cgs systems

        Returns:
            (np.ndarray or float) -> The converted value(s)
        """

        # Prevent modification to the origanal data
        if isinstance(data, np.ndarray):
            data = data.copy()

        if not isinstance(from_unit, UnitSystem):
            raise TypeError("Paramiter \"from_unit\" was not of type \"UnitSystem\".")
        if not isinstance(to_unit, UnitSystem):
            raise TypeError("Paramiter \"to_unit\" was not of type \"UnitSystem\".")

        # Do nothing if both units are the same
        if from_unit != to_unit:

            if from_unit == UnitSystem.h_less_comoving_GADGET:
                if to_unit == UnitSystem.physical:
                    if None in (h, h_scale_exponent, a, a_scale_exponent):
                        raise ValueError("Values were not provided for the following paramiters: \"h\", \"h_scale_exponent\", \"a\", \"a_scale_exponent\"")
                    data *= np.power(h, h_scale_exponent) * np.power(a, a_scale_exponent)

                elif to_unit == UnitSystem.cgs:
                    if cgs_conversion_factor is None:
                        raise ValueError("Values were not provided for the following paramiters: \"h\", \"h_scale_exponent\", \"a\", \"a_scale_exponent\", \"cgs_conversion_factor\"")
                    data = UnitSystem.convert_data(data, from_unit, UnitSystem.physical, h, h_scale_exponent, a, a_scale_exponent, cgs_conversion_factor) * cgs_conversion_factor
                    

                elif to_unit == UnitSystem.physical_comoving:
                    if None in (h, h_scale_exponent):
                        raise ValueError("Values were not provided for the following paramiters: \"h\", \"h_scale_exponent\", \"a\", \"a_scale_exponent\", \"cgs_conversion_factor\"")
                    data *= np.power(h, h_scale_exponent)


                elif to_unit == UnitSystem.cgs_comoving:
                    if cgs_conversion_factor is None:
                        raise ValueError("Values were not provided for the following paramiters: \"h\", \"h_scale_exponent\", \"a\", \"a_scale_exponent\", \"cgs_conversion_factor\"")
                    data = UnitSystem.convert_data(data, from_unit, UnitSystem.physical_comoving, h, h_scale_exponent, a, a_scale_exponent, cgs_conversion_factor) * cgs_conversion_factor

                else:
                    raise ValueError("Value of \"UnitSystem\" object was unexpected.")

            elif to_unit == UnitSystem.h_less_comoving_GADGET:
                if from_unit == UnitSystem.physical:
                    if None in (h, h_scale_exponent, a, a_scale_exponent):
                        raise ValueError("Values were not provided for the following paramiters: \"h\", \"h_scale_exponent\", \"a\", \"a_scale_exponent\"")
                    data /= np.power(h, h_scale_exponent) * np.power(a, a_scale_exponent)

                elif from_unit == UnitSystem.cgs:
                    if None in (h, h_scale_exponent, a, a_scale_exponent, cgs_conversion_factor):
                        raise ValueError("Values were not provided for the following paramiters: \"cgs_conversion_factor\"")
                    data /= np.power(h, h_scale_exponent) * np.power(a, a_scale_exponent) * cgs_conversion_factor

                elif from_unit == UnitSystem.physical_comoving:
                    if None in (h, h_scale_exponent):
                        raise ValueError("Values were not provided for the following paramiters: \"cgs_conversion_factor\"")
                    data /= np.power(h, h_scale_exponent)

                elif from_unit == UnitSystem.cgs_comoving:
                    if None in (h, h_scale_exponent, cgs_conversion_factor):
                        raise ValueError("Values were not provided for the following paramiters: \"cgs_conversion_factor\"")
                    data /= np.power(h, h_scale_exponent) * cgs_conversion_factor

            else:
                data = UnitSystem.convert_data(UnitSystem.convert_data(data, from_unit, UnitSystem.h_less_comoving_GADGET, h, h_scale_exponent, a, a_scale_exponent, cgs_conversion_factor), UnitSystem.h_less_comoving_GADGET, to_unit, h, h_scale_exponent, a, a_scale_exponent, cgs_conversion_factor)
    
        return data

    @staticmethod
    def convert_mass_to_solar(data: np.ndarray):
        """
        Converts mass value(s) to solar masses. This assumes input values are represented by the 'UnitSystem.cgs' system.

        Paramiters:
            (np.ndarray or float) data -> The value(s) to be converted

        Returns:
            (np.ndarray or float) -> The converted value(s)
        """

        return data / SimulationConstants.get_constants()["SOLAR_MASS"]

    @staticmethod
    def convert_mass_from_solar(data: np.ndarray):
        """
        Converts mass value(s) to 'UnitSystem.cgs'. This assumes input values are in solar masses.

        Paramiters:
            (np.ndarray or float) data -> The value(s) to be converted

        Returns:
            (np.ndarray or float) -> The converted value(s)
        """

        return data * SimulationConstants.get_constants()["SOLAR_MASS"]
