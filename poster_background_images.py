import numpy as np
from scipy.stats import norm
from matplotlib import pyplot as plt
from PIL import Image

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
    if i == 0: continue

    snapshot = ParticleReadConversion_EagleSnapshot(tag, simulations[i], model, relitive_data_root)

    r_200 = load_catalouge_field("Group_R_Crit200", "FOF", tag, simulations[i], SimulationModels.RECAL, relitive_data_root, UnitSystem.cgs)[assembily_history[simulations[i]][tag]["halo"] - 1]
    #selection_radius = r_200
    selection_radius = 9.257*(10**22)# 30 KPc

    particle_types = (ParticleType.gas, ParticleType.star, ParticleType.black_hole)
    net_angular_momentum = np.array([0, 0, 0], dtype = np.float64)
    for j in range(len(particle_types)):
        _, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_types[j], "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
        particle_velocities = snapshot.particle_read_sphere(particle_types[j], "Velocity", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
        #particle_velocities -= particle_velocities.mean(axis = 0)
        particle_masses = snapshot.particle_read_sphere(particle_types[j], "Mass", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)

        #net_angular_momentum += np.sum(np.cross(particle_locations_box_adjusted, particle_velocities) * particle_masses[:, None], axis = 0)
        net_angular_momentum += angular_momentum(particle_locations_box_adjusted, particle_velocities, particle_masses)

    rotation_axis_vector = net_angular_momentum / np.linalg.norm(net_angular_momentum, 2)
    rotation_plane_x_vector = produce_disk_x_basis(rotation_axis_vector)

    #particle_types = (ParticleType.gas, ParticleType.star)
    #particle_type_graph_names = ("Gas", "Star")
    #for j in range(len(particle_types)):
    #    particle_locations, particle_locations_box_adjusted = snapshot.particle_read_sphere(particle_types[j], "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    #    particle_locations_box_adjusted /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3

    #    particle_locations_transformed, _ = transform_coordinates(particle_locations_box_adjusted, rotation_plane_x_vector, rotation_axis_vector)

    #    fig = plt.figure(figsize = (24, 24))

    #    #ax = plt.subplot(2, 2, 3, aspect = "equal")
    #    ax = plt.gca()
    #    ax.set_facecolor("black")
    #    ax.set_xlabel("X (kpc)")
    #    ax.set_ylabel("Z (kpc)")
    #    ax.scatter(particle_locations_transformed[:, 0], particle_locations_transformed[:, 2], s = 0.1, c = "y")
    #    x_lims = ax.get_xlim()
    #    y_lims = ax.get_ylim()
    #    axis_lims = (min(x_lims[0], y_lims[0]), max(x_lims[1], y_lims[1]))
    #    ax.set_xlim(*axis_lims)
    #    ax.set_ylim(*axis_lims)
    
    #    fig.suptitle("{} Particle Projection of the Assembled Galaxy in the {} Regime".format(particle_type_graph_names[j], Simulations.to_string(simulations[i]).title()))
    #    plt.show()


    _, gas_particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.gas, "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    gas_particle_locations_box_adjusted /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3
    gas_particle_locations_transformed, _ = transform_coordinates(gas_particle_locations_box_adjusted, rotation_plane_x_vector, rotation_axis_vector)
    gas_density = snapshot.particle_read_sphere(ParticleType.gas, "Density", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    gas_smoothing_length = snapshot.particle_read_sphere(ParticleType.gas, "SmoothingLength", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    gas_smoothing_length /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3
    #gas_smoothing_length *= 0.1

    _, star_particle_locations_box_adjusted = snapshot.particle_read_sphere(ParticleType.star, "Coordinates", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs, return_coordinates = True)
    star_particle_locations_box_adjusted /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3
    star_particle_locations_transformed, _ = transform_coordinates(star_particle_locations_box_adjusted, rotation_plane_x_vector, rotation_axis_vector)
    star_smoothing_length = snapshot.particle_read_sphere(ParticleType.star, "SmoothingLength", galaxy_centres[i], selection_radius, UnitSystem.cgs, UnitSystem.cgs)
    star_smoothing_length /= constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3
    #star_smoothing_length *= 0.1

    #print(max(np.sqrt(np.sum(gas_particle_locations_transformed**2, axis = 1))))
    #print(max(gas_smoothing_length))
    #print()
    #print(max(np.sqrt(np.sum(star_particle_locations_transformed**2, axis = 1))))
    #print(max(star_smoothing_length))
    #print()
    #print()

    gas_xy = gas_particle_locations_transformed[:, 0:3:2]
    star_xy = star_particle_locations_transformed[:, 0:3:2]

    #aspect_ratio = 1
    aspect_ratio = 41 / 58

    width = 2000
    height = int(aspect_ratio * width)
    pixels = np.zeros(shape = (width, height), dtype = float)

    for x in range(width):
        #x_physical = (x - (width / 2)) * selection_radius / (constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3)
        x_physical = ((x / width) - 0.5) * 30
        for y in range(height):
            #y_physical = (y - (height / 2)) * selection_radius / (constants.SimulationConstants.get_constants()["CM_PER_MPC"] * 10**-3)
            y_physical = ((y / height) - 0.5) * 30 * aspect_ratio
            
            gas_intersection_filter = np.sqrt(np.sum((np.array([x_physical, y_physical]) - gas_xy)**2, axis = 1)) < gas_smoothing_length
            #gas_value = sum(gas_smoothing_length[gas_intersection_filter] - np.sqrt(np.sum((np.array([x_physical, y_physical]) - gas_xy[gas_intersection_filter])**2, axis = 1)))
            #gas_value = sum(np.sqrt(np.sum((np.array([x_physical, y_physical]) - gas_xy[gas_intersection_filter])**2, axis = 1)) / gas_smoothing_length[gas_intersection_filter])
            gas_value = sum(norm.pdf(np.sqrt(np.sum((np.array([x_physical, y_physical]) - gas_xy[gas_intersection_filter])**2, axis = 1)) / gas_smoothing_length[gas_intersection_filter], loc = 0, scale = 0.1))

            star_intersection_filter = np.sqrt(np.sum((np.array([x_physical, y_physical]) - star_xy)**2, axis = 1)) < star_smoothing_length
            #star_value = sum(star_smoothing_length[star_intersection_filter] - np.sqrt(np.sum((np.array([x_physical, y_physical]) - star_xy[star_intersection_filter])**2, axis = 1)))
            #star_value = sum(np.sqrt(np.sum((np.array([x_physical, y_physical]) - star_xy[star_intersection_filter])**2, axis = 1)) / star_smoothing_length[star_intersection_filter])
            star_value = sum(norm.pdf(np.sqrt(np.sum((np.array([x_physical, y_physical]) - star_xy[star_intersection_filter])**2, axis = 1)) / star_smoothing_length[star_intersection_filter], loc = 0, scale = 0.1))
            
            pixels[x, y] = star_value - gas_value
            if pixels[x, y] < 0:
                pixels[x, y] = 0
            
            if gas_value > 0 or star_value > 0: print("({}, {}) at ({:.2f}, {:.2f}):    {:.4f} - {:.4f}".format(x, y, x_physical, y_physical, star_value, gas_value))
                
    pixels *= 255 / pixels.max() * 0.7
    pixels = np.array(pixels, dtype = int)

    image = Image.fromarray(pixels.transpose()).convert("RGB")
    img_pixels = image.load()
    for pixel_x in range(image.size[0]):
        for pixel_y in range(image.size[1]):
            img_pixels[pixel_x, pixel_y] = (img_pixels[pixel_x, pixel_y][0], int(img_pixels[pixel_x, pixel_y][1] * 0.7), 0)
    image.save("{} rendered image.png".format(Simulations.to_string(simulations[i])))

    #print(pixels)
    #plt.imshow(pixels, cmap = "gray", vmin = 0, vmax = pixels.max())
    #plt.show()
