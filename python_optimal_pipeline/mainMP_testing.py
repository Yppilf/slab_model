import numpy as np
from convolver import Convolver
from molecule import Molecule
import os, copy
from multiprocessing import Pool
import pandas as pd # type: ignore

pd.set_option('future.no_silent_downcasting', True)

parameters_file = "data/broad_so2_dataset.txt"
convolver_file = "data/convolver_data/JWST_MIRI_MRS.json"
molecule_list = ["H2O", "OH", "CO", "CO2", "NH3", "HCN", "C2H2", "SO2"]
convolver_settings = ["lower", "upper", "optimal", "minimal"]
storagePath = "/scratch/s4950836/spectra2"
num_cores = 4

def genConvolvers(settings, convFile):
    convolvers = {}
    for setting in settings:
        conv = Convolver()
        conv.readResolvingModel(convFile, wavelength_overlap=setting)
        convolvers[setting] = conv
    return convolvers

def subprocess(params):
    slabs, storagePath, mol, temp, log10dens, convSetting, conv = params
    fname = f"{mol.molecule}_{temp}_{log10dens}_{convSetting}"
    fullPath = f"{storagePath}/{fname}.npz"

    if fullPath in slabs:
        return fullPath  # Slab already exists, no need to regenerate
    slabs.append(fullPath)  # Add new slab path

    newWavelengths = []
    newIntensities = []
    
    for j, channel in enumerate(conv.data):
        # For each wavelength range based on the convolver, generate a slab
        molec = copy.deepcopy(mol)
        molec.generateSlab(10**log10dens, temp, channel["wl"], channel["wu"])
        molec = conv.convolveData(molec, channel["wl"], channel["wu"])
        newWavelengths.extend(molec.convWavelength)
        newIntensities.extend(molec.convIntensity)

    # Sort wavelengths and intensities together
    sorted_indices = np.argsort(newWavelengths)
    newWavelengths = np.array(newWavelengths)[sorted_indices]
    newIntensities = np.array(newIntensities)[sorted_indices]

    # Save the slab data to a compressed file
    print(f"Writing to {fullPath}")
    np.savez_compressed(fullPath, wavelengths=newWavelengths, intensities=newIntensities)

    return fullPath

def genSlab(slabs, params, molecule_list, storagePath, convolvers):
    # Create a new set of parameters for (molecule, temperature, column_density, convolver_settings)
    conv = convolvers[params[-1]]

    # slabPaths = []
    # for i, mol in enumerate(molecule_list):
    #     fullPath = subprocess(slabs, mol, params[-2], params[i], params[-1], conv)
    #     slabPaths.append(fullPath)
    #     slabs.append(fullPath)

    param_list = [(slabs, storagePath, mol, params[-2], params[i], params[-1], conv) for i, mol in enumerate(molecule_list)]
    with Pool(num_cores) as pool:
        slabPaths = pool.map(subprocess, param_list)

    # slabPaths = [subprocess(slabs, storagePath, mol, params[-2], params[i], params[-1], conv) for i, mol in enumerate(molecule_list)]
    slabs.extend(slabPaths)

    return slabs, slabPaths
            
def combineSingleSlabs(params):
    filepaths, storagePath, idx = params
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

def main(convolver_settings, convolver_file, molecule_list, param_filepath, storagePath, num_cores):
    # Create a semaphore with a count of 1

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
    param_list = [(collection,storagePath,i) for i,collection in enumerate(singleMoleculeSlabs)]
    with Pool(num_cores) as pool:
        pool.map(combineSingleSlabs, param_list)

    # for i,collection in enumerate(singleMoleculeSlabs):
    #     combineSingleSlabs(collection, storagePath, i)

    # Phase 5
    for file in generated:
        os.remove(file)

if __name__ == "__main__":
    main(convolver_settings, convolver_file, molecule_list, parameters_file, storagePath, num_cores)