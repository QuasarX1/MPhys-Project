import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
from PIL import Image

import sys
sys.path.append("C:\\Scripts\\")
from cuda_factory import CUDA_Kernel

from DataAccess import ParticleReadConversion_EagleSnapshot, load_catalouge_field, Simulations, SimulationModels, ParticleType, UnitSystem, constants
from Physics import angular_momentum

from assembily_history import physical_centeral_mass_positions, assembily_history
from coordinate_transform import produce_disk_x_basis, transform_coordinates

model = SimulationModels.RECAL
tag = "028_z000p000"

relitive_data_root = ".\\gm_for_mphys"

galaxy_centres = [physical_centeral_mass_positions[sim][constants.tags[-1]] for sim in Simulations]#TODO: check whether cgs or physical code is uncommented!!!

simulations = (Simulations.Early, Simulations.Organic, Simulations.Late)
for i in range(len(simulations)):
    #if i == 0: continue

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulations[i], model, relitive_data_root)

    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulations[i], SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulations[i]][tag]["halo"] - 1]
    #selection_radius = r_200
    selection_radius = 9.257*(10**22)# 30 KPc

    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
    net_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    for j in range(len(particle_types)):
        _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_types[j], "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_velocities = snapshot.particle_read_sphere(particle_types[j], "Velocity", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        particle_masses = snapshot.particle_read_sphere(particle_types[j], "Mass", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        
        net_angular_momentum += angular_momentum(particle_locations_box_adjusted, particle_velocities, particle_masses)

    rotation_axis_vector = net_angular_momentum / np.linalg.norm(net_angular_momentum, 2)
    rotation_plane_x_vector = produce_disk_x_basis(rotation_axis_vector)



    _, gas_particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.gas, "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    gas_particle_locations_box_adjusted /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3
    gas_particle_locations_transformed, _ = transform_coordinates(gas_particle_locations_box_adjusted, rotation_plane_x_vector, rotation_axis_vector)
    gas_density = snapshot.particle_read_sphere(ParticleType.gas, "Density", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    gas_smoothing_length = snapshot.particle_read_sphere(ParticleType.gas, "SmoothingLength", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    gas_smoothing_length /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3

    _, star_particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.star, "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    star_particle_locations_box_adjusted /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3
    star_particle_locations_transformed, _ = transform_coordinates(star_particle_locations_box_adjusted, rotation_plane_x_vector, rotation_axis_vector)
    star_smoothing_length = snapshot.particle_read_sphere(ParticleType.star, "SmoothingLength", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    star_smoothing_length /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3



    #gas_xy = gas_particle_locations_transformed[:, 0:3:2]
    #star_xy = star_particle_locations_transformed[:, 0:3:2]
    gas_xy = gas_particle_locations_transformed[:, 1:3]
    star_xy = star_particle_locations_transformed[:, 1:3]

    aspect_ratio = 1
    #aspect_ratio = 41 / 58

    width = 2000
    height = int(aspect_ratio * width)
    
    c_func = """\
    extern "C" __global__
    void make_image(float *pixels, int width, float aspect_ratio, float physical_radius, int n_gas, int n_stars, float *gas_x, float *gas_y, float *star_x, float *star_y, float *gas_smoothing_length, float *star_smoothing_length)
    {
        int i;
        
        size_t tid = blockIdx.x * blockDim.x + threadIdx.x;

        if (tid == 0) { printf("%d, %f %d, %d\\n", width, aspect_ratio, n_gas, n_stars); }

        int height = (int)(aspect_ratio * width);

        if (tid < width * height)
        {
            int x = tid % width;
            float x_physical = (((float)x / (float)width) - 0.5) * physical_radius;

            int y = (int)(tid / width);
            float y_physical = (((float)y / (float)height) - 0.5) * physical_radius * aspect_ratio;

            // Gas
            float gas_value = 0.0;
            for (i = 0; i < n_gas; ++i)
            {
                float dist_to_particle = sqrtf(powf(x_physical - gas_x[i], 2.0) + powf(y_physical - gas_y[i], 2.0));
                if (dist_to_particle < gas_smoothing_length[i])
                {
                    float norm_x = dist_to_particle / gas_smoothing_length[i], norm_mu = 0.0, norm_sigma = 0.06, norm_a = 2.5;
                    //gas_value += 1.0;
                    gas_value += norm_a * powf(2.718281828459, -0.5 * powf((norm_x - norm_mu) / norm_sigma, 2.0)) / (norm_sigma * sqrtf(2.0 * 3.1415926535));
                }
            }

            // Star
            float star_value = 0.0;
            for (i = 0; i < n_stars; ++i)
            {
                float dist_to_particle = sqrtf(powf(x_physical - star_x[i], 2.0) + powf(y_physical - star_y[i], 2.0));
                if (dist_to_particle < star_smoothing_length[i])
                {
                    float norm_x = dist_to_particle / star_smoothing_length[i], norm_mu = 0.0, norm_sigma = 0.01, norm_a = 1.0;
                    //star_value += 1.0;
                    star_value += norm_a * powf(2.718281828459, -0.5 * powf((norm_x - norm_mu) / norm_sigma, 2.0)) / (norm_sigma * sqrtf(2.0 * 3.1415926535));
                }
            }

            pixels[tid] = star_value - gas_value;
            //pixels[tid] = star_value;
            //pixels[tid] = gas_value;
            if (pixels[tid] < 0.0) { pixels[tid] = 0.0; }
            printf("%f %f %f\\n", pixels[tid], star_value, gas_value);
        }
    }
    """
    n_pixels = width * height
    r = selection_radius / (constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3)
    make_image = CUDA_Kernel(c_func, "make_image",
                             [np.float32, np.int32, np.float32, np.float32, np.int32, np.int32, np.float32,  np.float32,  np.float32,   np.float32,   np.float32,  np.float32],
                             [n_pixels,   1,        1,          1,          1,        1,        len(gas_xy), len(gas_xy), len(star_xy), len(star_xy), len(gas_xy), len(star_xy)],
                             [0])
    make_image.blocks = CUDA_Kernel.minimum_needed_blocks(n_pixels)
    #make_image.threads_per_block = 512
    pixels = make_image(width, aspect_ratio, r, len(gas_xy), len(star_xy), gas_xy[:, 0], gas_xy[:, 1], star_xy[:, 0], star_xy[:, 1], gas_smoothing_length, star_smoothing_length)[:n_pixels]
    pixels = pixels.reshape((width, height))
    
    pixels = (pixels * 255 / pixels.max()).astype(int).transpose()
    
    pixels = pixels[:, :, np.newaxis]
    rgb_pixels = np.append(pixels, 0.7 * pixels, axis = 2)
    rgb_pixels = np.append(rgb_pixels, 0 * pixels, axis = 2)
    pixels = rgb_pixels.astype("uint8")
    image = Image.fromarray(pixels, "RGB")

    image.save("rendered image {}.png".format(Simulations.to_string(simulations[i])))
