from coniii import *
import numpy as np
from scipy.spatial.distance import cdist
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
import random
import math

def generate_assembly_solve(N, params_assembly_num, params_assembly_density):

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
    def create_rough_assembly(neuron_coords, num_points, mean, std_dev):
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

        # Convert arrays to tuples for comparison
        a_tuples = set(map(tuple, a))
        b_tuples = set(map(tuple, b))
        # Find the intersection
        intersection = a_tuples & b_tuples
        # Return the Szymkiewicz-Simpson coefficient
        return len(intersection) / min(len(a_tuples), len(b_tuples))

    def generate_assembly(neuron_array, assemblies, i, coeff_data,count, mean, dev):
        
        
        while True:
            params_distruibution_centre = np.random.uniform(0, mean, size=2)
            params_distribution_std_dev = np.random.uniform(dev-2, dev+2, size=2)
            assembly = create_rough_assembly(neuron_array, params_assembly_density, params_distruibution_centre, 1)
            
            if len(assembly) > 1:
                if i == 0:
#                     print("Number of neurons in assembly: ", len(assembly))
#                     print("Neurons in assembly:\n", assembly)
#                     print("----------------------------------------------------------------------------------------")
                    assemblies.append(assembly)
                    break
                else:
                    for j in range(i):
                        #if simpson_coefficient_2d(assembly, assemblies[j]) > 0.35:
                        #For better calcuation, set one assembly has 2 units and another one have 3. Also these 2 assemblies have one and only one identical element
                        if (len(set(tuple(x) for x in assembly.tolist()) & set(tuple(y) for y in assemblies[j].tolist())) != 1) or (len(assembly)==len(assemblies[j])):
                            break
                        else:
                            coeff_data.append((j, i, simpson_coefficient_2d(assembly, assemblies[j])))
                            count += 1
                    if count == i:
#                         print("Number of neurons in assembly: ", len(assembly))
#                         print("Neurons in assembly:\n", assembly)
#                         print("----------------------------------------------------------------------------------------")
#                         print("The %s and the %s assembly mean pairwise Szymkiewicz-Simpson coefficient value: %.2f"
#                               % (coeff_data[0][0], coeff_data[0][1], float(coeff_data[0][2])))
#                         print("Which is less than 0.35")
                        assemblies.append(assembly)
                        count = 0
                        break
    
    bpmean_para = (round(math.sqrt(N))  - 1)
    bpdev_para = round((round(math.sqrt(N)))/2.0)
    
    count = 0
    neuron_array = create_neuron_array(N)
    assemblies = []
    coeff_data = []
    for i in range(params_assembly_num):
        generate_assembly(neuron_array, assemblies, i, coeff_data, count, bpmean_para, bpdev_para)
        
    coord_to_index = {tuple(coord): i for i, coord in enumerate(neuron_array)}
    neuron_array_coord = neuron_array
    assemblies_coord = assemblies
    neuron_array = np.array([coord_to_index[tuple(coord)] for coord in neuron_array])
    for i in range(len(assemblies)):
        assemblies[i] = [coord_to_index[tuple(coord)] for coord in assemblies_coord[i] if tuple(coord) in coord_to_index]
    
    return assemblies