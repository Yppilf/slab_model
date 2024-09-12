import numpy as np
from slab_model import generateSlabModel
import pandas as pd
pd.set_option('future.no_silent_downcasting', True)

parameters_file = "permutations3.txt"

molecule_list = ["H2O", "OH", "CO", "CO2", "NH3", "SO2"] 

def save_spectrum(filename, wavelengths, summed_intensities, numberdensities):
    # Save the wavelengths, summed intensities, and number densities as compressed .npy
    np.savez_compressed(filename, wavelengths=wavelengths, intensities=summed_intensities, numberdensities=numberdensities)

def process_spectra(parameters_file, output_dir):
    with open(parameters_file, 'r') as file:
        for line_num, line in enumerate(file, start=1):
            if line_num == 1: 
                continue    # first spectrum is generated in testing areadt
            params = eval(line.strip())  # Converts the string of parameters back into a tuple/list
            print(f"Processing spectrum {line_num}")

            # Generate the spectrum for these parameters
            spectrum_list = generateSlabModel(molecule_list, params[-2], params[:-2], overlap=params[-1])

            # Assuming all elements in the spectrum have the same 'wavelengths'
            wavelengths = spectrum_list[0]["wavelengths"]

            # Sum the intensities for all molecules
            summed_intensities = np.zeros_like(wavelengths)
            numberdensities = {}

            for element in spectrum_list:
                summed_intensities += element["intensities"]
                numberdensities[element["molecule"]] = element["numberdensity"]

            # Save the wavelengths, summed intensities, and number densities to disk
            filename = f"{output_dir}/spectrum_{line_num}.npz"
            save_spectrum(filename, wavelengths, summed_intensities, numberdensities)

            break

process_spectra(parameters_file, "/scratch/s4950836/spectra")