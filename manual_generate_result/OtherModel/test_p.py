import os
os.environ['MPLCONFIGDIR'] = '/tmp/mpl_config'

import warnings
# 仅在这个context里忽略警告
with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    from coniii import *
import numpy as np
from tqdm.auto import tqdm
from scipy.spatial.distance import cdist
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt
import random
import Jangenerate_assembly #In the same dir
import Jangenerate_SpikeCount #In the same dir
from scipy.stats import poisson

import itertools
import time
import math

warnings.filterwarnings("ignore")

# Define Parameters
T = 3600 # time of simul“ation
dT = 0.5 # time step
params_assembly_num =5 # number of assemblies
params_point_into_neuron_distance = 0.5 

# Length of an active event as a number of timesteps
eventDur = np.random.randint(1, 10)
# Probability with which a unit is particularly active in a single timestep
eventProb = np.random.uniform(0.01, 0.05)
# Firing rate multiplier at active events
eventMult = np.random.uniform(6, 10)  # random number between 1 and 5
showPlot = False

def binaryOutput(original_list):
    # Create a new list to hold the tuples
    tuples_list = []

    # Generate all possible combinations of two elements for each sublist
    for sublist in original_list:
        combinations = itertools.combinations(sublist, 2)
        # Convert the combinations into tuples and add them to the list
        tuples_list.extend(tuple(sorted(combination)) for combination in combinations)

    # Remove duplicates by converting the list to a set then back to a list
    unique_tuples = list(set(tuples_list))
    
    return unique_tuples

N = [25,36]
params_assembly_density = [n // 6 for n in N] # size of neurons in each assembly

assemblies_list = []
spikeCount_list = []
binary_list = []

for i in tqdm(range (0, len(N))):
    print("1")
    fire_rate_background = np.random.uniform(1, 6, N[i])
    assemblies = Jangenerate_assembly.generate_assembly_solve(N[i], params_assembly_num, params_assembly_density[i])
    # Output 0, 1 type spikes
    spikeCount = Jangenerate_SpikeCount.generateSpikeCountSolve(N[i], T, dT, assemblies, (1, 6), eventDur, eventProb, eventMult, showPlot)
    # Transform to -1, 1 distribution
    spikeCount[spikeCount == 0] = -1
    assemblies_list.append(assemblies)
    spikeCount_list.append(spikeCount)
    print(assemblies)
    print("_______________________________________")
    binary_list.append(binaryOutput(assemblies))