import h5py as h5

tags = ("000_z020p000", "001_z015p132", "002_z009p993",
        "003_z008p988", "004_z008p075", "005_z007p050",
        "006_z005p971", "007_z005p487", "008_z005p037",
        "009_z004p485", "010_z003p984", "011_z003p528",
        "012_z003p017", "013_z002p478", "014_z002p237",
        "015_z002p012", "016_z001p737", "017_z001p487",
        "018_z001p259", "019_z001p004", "020_z000p865",
        "021_z000p736", "022_z000p615", "023_z000p503",
        "024_z000p366", "025_z000p271", "026_z000p183",
        "027_z000p101", "028_z000p000")

class SimulationConstants(object):
    __simulation_constants = None

    def __init__(self):
        #from DataAccess._catalouge_access import load_catalouge_field
        #load_catalouge_field
        
        with h5.File("./gm_for_mphys/RECAL/Organic/groups_028_z000p000/eagle_subfind_tab_028_z000p000.0.hdf5", "r") as datafile:
            self.__values = dict(datafile["Constants"].attrs)

    @staticmethod
    def get_constants():
        if SimulationConstants.__simulation_constants is None:
            SimulationConstants.__simulation_constants = SimulationConstants()
            
        return SimulationConstants.__simulation_constants.__values.copy()
