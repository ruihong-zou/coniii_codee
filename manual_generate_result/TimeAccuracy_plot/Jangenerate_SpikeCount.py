import numpy as np
from scipy.stats import poisson
import matplotlib.pyplot as plt

def generateSpikeCountSolve(N, T, dT, assemblies, rateRange, eventDur, eventProb, eventMult, showPlot):
    
    """
    Function to generate spike count for a given set of parameters.
    N: Number of neurons
    T: Number of time steps
    dT: Width of a time step in units of seconds
    assemblies: List of collections of units to be combined to assemblies
    rateRange: Range of base firing rate in units of inverse seconds (Hertz)
    eventDur: Length of an active event as a number of timesteps
    eventProb: Probability with which a unit is particularly active in a single timestep
    eventMult: Firing rate multiplier at active events
    """
    
    # Initialization
    fire_rate_background = np.random.uniform(rateRange[0], rateRange[1], size=(T, N))
    activation_field = np.zeros((T, N), dtype=bool)

    # Generate activation for each assembly and each neuron
    for t in range(T):
        for assembly in assemblies:
            # If any neuron in this assembly is activated, all neurons in this assembly are activated 
            if np.any(np.random.rand(len(assembly)) < eventProb):
                activation_field[t, assembly] = True
        # For neurons not in any assembly, they are activated independently
        not_in_assembly = np.setdiff1d(np.arange(N), np.concatenate(assemblies))
        activation_field[t, not_in_assembly] = np.random.rand(len(not_in_assembly)) < eventProb

    # If neuron is activated at a timestep, its firing rate is increased by eventMult times
    fire_rate = np.where(activation_field, fire_rate_background * eventMult, fire_rate_background)
    # Generate Poisson spike count for each neuron at each timestep
    spike_count = np.random.poisson(fire_rate * dT)
    spike_binary = np.where(spike_count < 6, 0, 1)
    
    if showPlot:
        plt.figure(figsize=(8, 6), dpi=80)
        plt.imshow(spike_binary.T, aspect="auto", cmap="gray_r", interpolation="none")
    
    return spike_binary