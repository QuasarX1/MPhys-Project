"""
Provides functional and OOP access to the EAGLE GM simulation data.

Code is based on the tutorial by J. Davies avalible at "https://j-davies-ari.github.io/eagle-guide/".
"""

import DataAccess._constants as constants
from DataAccess._unit_conversions import UnitSystem
from DataAccess._particle_type import ParticleType, all_particle_types, baryonic_particle_types, baryonic_low_mass_particle_types
from DataAccess._simulation_combinations import Simulations, SimulationModels
from DataAccess._catalouge_access import load_catalouge_field, particles_by_group
from DataAccess._snapshots import ParticleReadConversion_EagleSnapshot
