from enum import Enum

class Simulations(Enum):
    """
    Simulations that are accepted paramiters for functions in the "DataAccess" package.
    """
    Early = 0
    Organic = 1
    Late = 2

    @staticmethod
    def to_string(value) -> str:
        """
        Returns the appropriate string representation of the "Simulations" enum's values.

        Paramiters:
            Simulations value -> The enum value to convert

        Returns:
            str -> The coresponding string representation
        """

        if not isinstance(value, Simulations):
            raise TypeError("Paramiter \"value\" was not of type \"Simulations\".")

        if value == Simulations.Early:
            return "Early"
        elif value == Simulations.Organic:
            return "Organic"
        elif value == Simulations.Late:
            return "Late"
        else:
            raise ValueError("Value of \"Simulations\" object was unexpected.")

    @staticmethod
    def from_string(value: str):
        """
        Returns the appropriate "Simulations" enum value from the coresponding string representation.

        Paramiters:
            str value -> The string to be converted

        Returns:
            Simulations -> The coresponding enum value
        """

        if value not in ("Early", "Organic", "Late"):
            raise ValueError("Paramiter \"value\" was not of type \"Simulations\".")

        if value == "Early":
            return Simulations.Early
        elif value == "Organic":
            return Simulations.Organic
        elif value == "Late":
            return Simulations.Late
        else:
            raise ValueError("Value of string representation was unexpected.")



class SimulationModels(Enum):
    """
    Simulation models that are accepted paramiters for functions in the "DataAccess" package.
    """
    RECAL = 0
    DMONLY = 1

    @staticmethod
    def to_string(value) -> str:
        """
        Returns the appropriate string representation of the "SimulationModels" enum's values.

        Paramiters:
            SimulationModels value -> The enum value to convert

        Returns:
            str -> The coresponding string representation
        """

        if not isinstance(value, SimulationModels):
            raise TypeError("Paramiter \"value\" was not of type \"SimulationModels\".")

        if value == SimulationModels.RECAL:
            return "RECAL"
        elif value == SimulationModels.DMONLY:
            return "DMONLY"
        else:
            raise ValueError("Value of \"SimulationModels\" object was unexpected.")

    @staticmethod
    def from_string(value: str):
        """
        Returns the appropriate "SimulationModels" enum value from the coresponding string representation.

        Paramiters:
            str value -> The string to be converted

        Returns:
            SimulationModels -> The coresponding enum value
        """

        if value not in ("RECAL", "DMONLY"):
            raise ValueError("Paramiter \"value\" was not of type \"SimulationModels\".")

        if value == "RECAL":
            return SimulationModels.RECAL
        elif value == "DMONLY":
            return SimulationModels.DMONLY
        else:
            raise ValueError("Value of string representation was unexpected.")
