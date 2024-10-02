import numpy as np
import random

minTemp = 100   # Minimum temperature of the slab model in Kelvin
maxTemp = 1200  # Maximum temperature of the slab model in Kelvin
dTemp = 100

molecule_list = ["H2O", "OH", "CO", "CO2", "NH3", "HCN", "C2H2"]                    # List of molecules to include in the data (Only most abundant species )
log10_densities = [18,15,18,14,17,14,14]                                            # Log10 of estimate of common column densities in cm^-2
molecule_target = "SO2"
target_min_abundance = 1
target_max_abundance = 18
dtarget_abundance = 1

convolver_settings = ["lower", "upper", "optimal", "minimal"]

temp_range = np.arange(minTemp, maxTemp, dTemp)
density_range = np.arange(target_min_abundance, target_max_abundance, dtarget_abundance)

paramList = []
for temp in temp_range:
    for dens in density_range:
        convSetting = random.choice(convolver_settings)
        params = (*log10_densities, dens, temp, convSetting)
        paramList.append(params)

with open('permutations4.txt', 'w') as file:
    for combi in paramList:
        file.write(f"{combi}\n")
