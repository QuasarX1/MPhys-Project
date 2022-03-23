from enum import Enum

class ParticleType(Enum):
    """
    GADGET particle types that are accepted paramiters for functions in the "DataAccess" package.

    Use ParticleType.<option>.value to get the integer associated with the particle type.
    """
    gas = 0
    dark_matter = 1
    star = 4
    black_hole = 5

    @staticmethod
    def to_string(value) -> str:
        """
        Returns the appropriate string representation of the "ParticleType" enum's values.

        Paramiters:
            ParticleType value -> The enum value to convert

        Returns:
            str -> The coresponding string representation
        """

        if not isinstance(value, ParticleType):
            raise TypeError("Paramiter \"value\" was not of type \"ParticleType\".")

        if value == ParticleType.gas:
            return "Gas"
        elif value == ParticleType.dark_matter:
            return "Dark Matter"
        elif value == ParticleType.star:
            return "Star"
        elif value == ParticleType.black_hole:
            return "Black Hole"
        else:
            raise ValueError("Value of \"ParticleType\" object was unexpected.")

    @staticmethod
    def from_string(value: str):
        """
        Returns the appropriate "ParticleType" enum value from the coresponding string representation.

        Paramiters:
            str value -> The string to be converted

        Returns:
            ParticleType -> The coresponding enum value
        """

        if value not in ("Gas", "Dark Matter", "Star", "Black Hole"):
            raise ValueError("Paramiter \"value\" was not of type \"Simulations\".")

        if value == "Gas":
            return ParticleType.gas
        elif value == "Dark Matter":
            return ParticleType.dark_matter
        elif value == "Star":
            return ParticleType.star
        elif value == "Black Hole":
            return ParticleType.black_hole
        else:
            raise ValueError("Value of string representation was unexpected.")

all_particle_types = (ParticleType.gas, ParticleType.dark_matter, ParticleType.star, ParticleType.black_hole)
barionic_particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
barionic_low_mass_particle_types = (ParticleType.gas, ParticleType.star)
