import numpy as np
import itertools, random

### User definitions

minTemp = 100   # Minimum temperature of the slab model in Kelvin
maxTemp = 2000  # Maximum temperature of the slab model in Kelvin
flipTemp = 600  # Temperature from which point we will use a different temperature step size
dTemp1 = 200    # Temperature step size in Kelvin before flip   
dTemp2 = 100    # Temperature step size in Kelvin after flip

molecule_list = ["H2O", "OH", "CO", "CO2", "NH3", "SO2"]                    # List of molecules to include in the data (Only most abundant species )
log10_lowbound_densities = [18,15,18,14,17,14]                              # Log10 of lower estimate of common column densities in cm^-2
log10_highbound_densities = [18,15,19,16,18,15]                             # Log10 of upper estimate of common column densities in cm^-2
density_focus = [1, 1, 2, 2, 1, 0.5]                                        # Step sizes taken in log10 space in column density, allowing for focussing on specific species and less so on less important species
density_padding = 0                                                         # Extra padding around common ranges in column density in log10
low_presence_range = [1,2]                                                  # Additional range for low presence of molecules in log10
low_density_treshold = 10                                                   # Remove spectra with multiple elements having a density below this treshold

### End of User definitions

### Program definitions 

convolver_settings = ["lower", "upper", "optimal", "minimal"]

### End of program definitions

# Calculate lists of parameters
tempRange = np.concatenate((
    np.arange(minTemp, flipTemp, dTemp1, dtype=int),    
    np.arange(flipTemp, maxTemp, dTemp2, dtype=int)
))
densityLists = []

for m,l,u,dd in zip(molecule_list,log10_lowbound_densities,log10_highbound_densities,density_focus):
    densityList = np.arange(l-density_padding, u+density_padding+dd, dd, dtype=float)
    densityList = np.unique(np.concatenate((low_presence_range, densityList)))
    densityLists.append(densityList)

all_lists = densityLists + [tempRange]

total_combinations = list(itertools.product(*all_lists))

# Filter based on the low-density threshold
filtered_combinations = []
for combination in total_combinations:
    low_density_count = sum([density < low_density_treshold for density in combination[:len(molecule_list)]])
    if low_density_count <= 1:
        filtered_combinations.append(combination)

# Randomly assign one convolver setting per combination to avoid exponential growth
convolver_random = [random.choice(convolver_settings) for _ in filtered_combinations]
filtered_combinations = [(*comb, conv) for comb, conv in zip(filtered_combinations, convolver_random)]

with open('permutations3.txt', 'w') as file:
    for combi in filtered_combinations:
        file.write(f"{combi}\n")




