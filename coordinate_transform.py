import numpy as np

def transform_coordinates(coordinates: np.ndarray, e1: np.ndarray, e2: np.ndarray, e3: np.ndarray = None):
    """
    Take coordinates and re-express them based on a caretesian system with the same origin.

    p'_i = sum_j(p_j * e_j . e'_i)
    """

    provide_third_basis = False
    if e3 is None:
        provide_third_basis = True
        length_product = np.linalg.norm(e1, 2) * np.linalg.norm(e2, 2)
        e3 = np.cross(e1, e2) / (length_product * np.sin(np.arccos(e1.dot(e2) / length_product)))

    new_bases = (e1, e2, e3)

    new_coordinates = np.empty(coordinates.shape)

    for new_coordinate_index in range(len(new_coordinates)):
        for new_ordinate_index in range(len(new_coordinates[0])):
            new_coordinates[new_coordinate_index, new_ordinate_index] = 0
            for old_ordinate_index in range(len(new_coordinates[0])):
                new_coordinates[new_coordinate_index, new_ordinate_index] += coordinates[new_coordinate_index, old_ordinate_index] * new_bases[new_ordinate_index][old_ordinate_index]

    if provide_third_basis:
        return new_coordinates, e3
    else:
        return new_coordinates

def produce_disk_x_basis(axis_vector: np.ndarray):
    """
    Create an arbitrary primary coordinate basis given the axis of rotation in global cartesian coordinates.

    The returned vector is arbitraraly one of two possible vectors that lie purely in the global z=0 plane.
    """
    return np.array([axis_vector[1], -axis_vector[0], 0]) / np.sqrt(axis_vector[0]**2 + axis_vector[1]**2)
    #if axis_vector[1] != 0:
    #    projected_axis_ratio = axis_vector[0] / axis_vector[1]
    #    return np.array([1, -projected_axis_ratio, 0]) / np.sqrt(1 + (projected_axis_ratio)**2)
    #else:
    #    return np.array([0, -1, 0])



if __name__ == "__main__":
    # Define new bases in terms of current cartesian coords
    e1 = np.array([0, -1, 0])
    e2 = np.array([-1, 0, 0])
    e3 = np.array([0, 0, -1])

    # Define coordinates in current cartesian coords
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])

    # Convert
    print(transform_coordinates(coords, e1, e2, e3))



    from matplotlib import pyplot as plt

    e1 = np.array([2**-0.5, 2**-0.5, 0])
    e2 = np.array([-(2**-0.5), 2**-0.5, 0])
    e3 = np.array([0, 0, 1])

    #e1 = np.array([3**-0.5, 3**-0.5, -(3**-0.5)])
    #e2 = np.array([-(3**-0.5), 3**-0.5, -(3**-0.5)])

    pattern_x = []
    pattern_y = []
    n = 100
    for i in range(n):
        pattern_x.append(np.random.randint(-1000, 1000))
        pattern_y.append(np.random.randint(-1000, 1000))

    #pattern_x = [-1, 0, 1, 1, 0]
    #pattern_y = [-1, 0, 0, 1, 1]
    #n = len(pattern_x)

    new_coords, new_basis = transform_coordinates(np.array([[pattern_x[i], pattern_y[i], 0] for i in range(n)]), e1, e2)
    print("New basis =", new_basis, "expected:",e3 )

    new_x = [new_coords[i][0] for i in range(n)]
    new_y = [new_coords[i][1] for i in range(n)]

    fig = plt.figure()
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.scatter(pattern_x, pattern_y)
    ax1.set_aspect('equal', adjustable='box')
    ax2 = fig.add_subplot(1, 2, 2)
    ax2.scatter(new_x, new_y)
    ax2.set_aspect('equal', adjustable='box')

    lower_x = min(ax1.get_xlim()[0], ax2.get_xlim()[0])
    upper_x = max(ax1.get_xlim()[1], ax2.get_xlim()[1])
    lower_y = min(ax1.get_ylim()[0], ax2.get_ylim()[0])
    upper_y = max(ax1.get_ylim()[1], ax2.get_ylim()[1])
    ax1.set_xlim(lower_x, upper_x)
    ax1.set_ylim(lower_y, upper_y)
    ax2.set_xlim(lower_x, upper_x)
    ax2.set_ylim(lower_y, upper_y)
    plt.show()


    #import time
    #start = time.time()
    #for _ in range(10000000):
    #    produce_disk_x_basis(np.array([2**-0.5, 2**-0.5, 0]))
    #print(time.time() - start)

    print(produce_disk_x_basis(np.array([0, 1, 0])))
    print(produce_disk_x_basis(np.array([1, 0, 0])))
    print(produce_disk_x_basis(np.array([2**-0.5, 2**-0.5, 0])))
    print(produce_disk_x_basis(np.array([-(2**-0.5), 2**-0.5, 0])))