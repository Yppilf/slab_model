import numpy as np
from convolver import Convolver
from molecule import Molecule
import os, copy
import pandas as pd # type: ignore
pd.set_option('future.no_silent_downcasting', True)

parameters_file = "data/testPermutations.txt"
convolver_file = "data/convolver_data/JWST_MIRI_MRS.json"
molecule_list = ["H2O", "OH", "CO", "CO2", "NH3", "SO2"] 
convolver_settings = ["lower", "upper", "optimal", "minimal"]
wavelengthRange = (5.00, 27.50)
# storagePath = "/scratch/s4950836"
storagePath = "."

def genConvolvers(settings, convFile):
    convolvers = {}
    for setting in settings:
        conv = Convolver()
        conv.readResolvingModel(convFile, wavelength_overlap=setting)
        convolvers[setting] = conv
    return convolvers

def genSlab(slabs, params, molecule_list, storagePath, convolvers):
    slabPaths = []

    # Create a new set of parameters for (molecule, temperature, column_density, convolver_settings)
    conv = convolvers[params[-1]]
    for i,mol in enumerate(molecule_list):
        fname = f"{mol.molecule}_{params[-2]}_{params[i]}_{params[-1]}"
        fullPath = f"{storagePath}/{fname}.npz"
        slabPaths.append(fullPath)
        if fname in slabs:
            continue
        slabs.append(fullPath)

        newWavelenghts = []
        newIntensities = []
        for j,channel in enumerate(conv.data):
            # For each wavelengthrange /* 1 channel */ based on the convolver, generate a slab
            molec = copy.deepcopy(mol)
            molec.generateSlab(10**params[i], params[-2], channel["wl"], channel["wu"])
            molec = conv.convolveData(molec, channel["wl"], channel["wu"])
            newWavelenghts.extend(molec.convWavelength)
            newIntensities.extend(molec.convIntensity)

        np.savez_compressed(fullPath, wavelengths=newWavelenghts, intensities=newIntensities)

    return slabs, slabPaths
            
def combineSingleSlabs(filepaths, storagePath, idx):
    totalIntensity = None
    wavelengths = []
    molecules = []
    numberDensitites = []

    # Add the intensities together for all the molecules (Wavelengths are the same since same convolver is used)
    for path in filepaths:
        loadedSpectrum = np.load(path)
        params = path.split("/")[-1].split(".npz")[0].split("_")
        molecules.append(params[0])
        numberDensitites.append(float(params[2]))
        if len(wavelengths) < 1:
            wavelengths = loadedSpectrum["wavelengths"]
            totalIntensity = np.array(loadedSpectrum["intensities"])
        else:
            totalIntensity += loadedSpectrum["intensities"]


    # Normalize the spectrum to 1 (Not relevant for data generation but makes it easier later on to process)
    totalIntensityNormalized = totalIntensity/np.max(totalIntensity)

    np.savez_compressed(f"{storagePath}/spectrum{idx}.npz", wavelengths=wavelengths, intensities=totalIntensityNormalized, molecules=molecules, numberDensitites=numberDensitites)

def main(convolver_settings, convolver_file, molecule_list, param_filepath, wavelengthRange, storagePath):
    # Phase 1
    # Generate a convolver object for each type of setting present/possible
    convolvers = genConvolvers(convolver_settings, convolver_file)

    # For each molecule modelled, calculate the molar weight
    molecules = [Molecule(mol) for mol in molecule_list]

    # Read the parameters for which the slab models need to be generated, each line being one set of parameters
    singleMoleculeSlabs = []
    with open(param_filepath, 'r') as file:
        generated = []
        for line_num, line in enumerate(file, start=1):
            params = eval(line.strip())  # Converts the string of parameters back into a tuple/list
            
            # Phase 2,3
            generated, slabPaths = genSlab(generated, params, molecules, storagePath, convolvers)
            singleMoleculeSlabs.append(slabPaths)
    
    # Phase 4
    for i,collection in enumerate(singleMoleculeSlabs):
        combineSingleSlabs(collection, storagePath, i)

    # Phase 5
    # for file in generated:
    #     os.remove(file)

if __name__ == "__main__":
    main(convolver_settings, convolver_file, molecule_list, parameters_file, wavelengthRange, storagePath)