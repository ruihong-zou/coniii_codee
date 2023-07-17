from coniii import *
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
import random

# Define Parameters
T = 3600 # time of simulation
N = 25  # number of neurons
dT = 0.5 # time step
params_assembly_density = 5 # size of neurons in each assembly
params_assembly_num = 3 # number of assemblies
params_point_into_neuron_distance = 0.5
fire_rate_background = np.random.uniform(1, 6, N)

# Function to create a neuron array using unit-spacing hexagonal lattice points neuron
def create_neuron_array(num_neurons):
    coords = []
    num_rows = int(np.sqrt(num_neurons))
    num_cols = num_rows
    for i in range(num_rows):
        for j in range(num_cols):
            x = j + 0.5 * (i % 2)  # Stagger every other row to create hexagonal pattern
            y = i * np.sqrt(3) / 2  # Scale row spacing by sqrt(3)/2 to create hexagonal pattern
            coords.append((x, y))
    return np.array(coords)

# Function to create an assembly of neurons
def create_assembly(neuron_coords, num_points, mean, std_dev):
    # Create the covariance matrix
    cov_matrix = np.eye(2) * std_dev**2
    # Draw points from a two-dimensional normal distribution
    points = np.random.multivariate_normal(mean, cov_matrix, num_points)
    # Find the closest neuron to each drawn point
    distances = cdist(points, neuron_coords)
    closest_neurons = np.argmin(distances, axis=1)
    # If the distance from the neuron to the drawn point is less than 0.5, it is considered part of the assembly
    #Change the parameter because of small size of neurons
    assembly_neurons = neuron_coords[closest_neurons[distances[np.arange(len(closest_neurons)), closest_neurons] < 2]]
    new_assembly_neurons = np.unique(assembly_neurons, axis = 0)

    return new_assembly_neurons

def simpson_coefficient_2d(a, b):
    print(a)
    print(b)
    print(6666666666666666666666666)
    # Convert arrays to tuples for comparison
    a_tuples = set(map(tuple, a))
    b_tuples = set(map(tuple, b))
    # Find the intersection
    intersection = a_tuples & b_tuples
    # Return the Szymkiewicz-Simpson coefficient
    return len(intersection) / min(len(a_tuples), len(b_tuples))

# Generate the neuron array
neuron_array = create_neuron_array(N)
first_time = True

# Generate assemblies
assemblies = []
count = 0
coeff_data = []
for i in range(params_assembly_num):
    while True:
        # Choose a centre from a uniform distribution
        params_distruibution_centre = np.random.uniform(0, 4, size=2)
        # Choose a standard deviation from a uniform distribution
        params_distribution_std_dev = np.random.uniform(1, 3, size=2)
        assembly = create_assembly(neuron_array,params_assembly_density ,params_distruibution_centre, 1)
        # If assembly is not empty, append it to assemblies and break the while loop
        if len(assembly) > 1:
            if first_time:
                first_time = False
                print("Number of neurons in assembly: ", len(assembly))
                print("Neurons in assembly:\n ", assembly)
                assemblies.append(assembly)
                break
            else:
                for j in range(i):
                    if (simpson_coefficient_2d(assembly, assemblies[j]) > 0.35 ):
                        break
                    else:
                        coeff_data.append((j, i, simpson_coefficient_2d(assembly, assemblies[j])))
                        count += 1
                if(count == i):
                    print(coeff_data)
                    print("Number of neurons in assembly: ", len(assembly))
                    print("Neurons in assembly:\n ", assembly)
                    assemblies.append(assembly)
                    break


# Plot the neuron array and assemblies
plt.figure(figsize=(20, 20))
plt.scatter(neuron_array[:, 0], neuron_array[:, 1], color='black', s=10)
for assembly in assemblies:
    plt.scatter(assembly[:, 0], assembly[:, 1], color='red', s=10)
plt.show()