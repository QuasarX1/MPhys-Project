import numpy as np

def time_from_redshift(z, omega0, omega_lambda, h):
    if omega_lambda != 0 and np.abs(omega0 - 1.0 + omega_lambda) > 10**-6:
        print(f"Omega0 = {omega0}")
        print("Omega_tot = {}".format(omega0 + omega_lambda))
        print("Descrepency = {}".format(omega0 - 1 + omega_lambda))
        raise ValueError("Invalid inputs - see output for details.")

    hTime = 9.778 / h# Gyrs

    if omega_lambda == 0:
        if omega0 == 1:
            time = (2 / 3) * hTime / ((1 + z)**1.5)
        else:
            if omega0 > 1:
                omega0_z = omega0 * z
                factor1 = omega0 / (2 * (omega0 - 1)**1.5)
                term1 = np.arccos((omega0_z - omega0 + 2) / (omega0_z + omega0))
                term2 = 2 * np.sqrt((omega0 - 1) * (omega0_z + 1) / (omega0_z + omega0))
                factor2 = term1 - term2
                time = factor1 * factor2 * hTime

            else:
                raise ValueError("Failed tryint to intergrate for open, matter dominated cosmology.")

    else:
        factor1 = (2 / 3) / np.sqrt(omega_lambda)
        term1 = np.sqrt(omega_lambda / omega0) / ((1 + z)**1.5)
        term2 = np.sqrt(1 + (omega_lambda / (omega0 * (1 + z)**3)))
        factor2 = np.log(term1 + term2)
        time = factor1 * factor2 * hTime

    return time
