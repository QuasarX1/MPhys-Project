"""
Provides functional and OOP access to the EAGLE GM simulation data.

Code is based on the tutorial by J. Davies avalible at "https://j-davies-ari.github.io/eagle-guide/".
"""

from ._unit_conversions import UnitSystem
from ._particle_type import ParticleType
from ._simulation_combinations import Simulations, SimulationModels
from ._catalouge_access import load_catalouge_field
from ._snapshots import ParticleReadConversion_EagleSnapshot
