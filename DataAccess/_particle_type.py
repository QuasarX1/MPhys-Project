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
