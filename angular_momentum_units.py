from DataAccess import constants, UnitSystem

def angular_momentum_units(value):
    """
    g cm^2 s^-1 -> M_sun kPc km s^-1
    """
    return UnitSystem.convert_mass_to_solar(value) * 10**(-2) / constants.SimulationConstants.get_constants()["CM_PER_MPC"]

def specific_angular_momentum_units(value):
    """
    cm^2 s^-1 -> kPc km s^-1
    """
    return value * 10**(-2) / constants.SimulationConstants.get_constants()["CM_PER_MPC"]
