from matplotlib.pyplot import figure,show,savefig
import copy, json

from molecule import Molecule
from convolver import Convolver

def generateSlabModel(mols, temperature, numberDensities, verbose=False, saveFig=False, overlap="optimal"):
    # Generate convolver for wavelength ranges
    convolver = Convolver()
    convolver.readResolvingModel("./convolver_data/JWST_MIRI_MRS.json", wavelength_overlap=overlap)     # Maybe for the report, plot a mosaic of different convolver settings?

    # generate molecules
    molecules = [Molecule(m, temperature) for m in mols]

    # Generate data for each of the molecules for each of the channels
    spectrum = []
    for idx, m in enumerate(molecules):
        newWavelengths = []
        newIntensities = []
        for i,channel in enumerate(convolver.data):
            if verbose:
                print(f"{m.molecule}: Channel {i+1}/{len(convolver.data)}")
            molecule = copy.deepcopy(m)     # Not sure if this is necessary, but as to not overwrite the original m each time
            molecule.get_moldata(channel["wl"], channel["wu"])
            molecule.generateSlab(numberDensities[idx], channel["wl"], channel["wu"])
            molecule = convolver.convolveData(molecule, channel["wl"], channel["wu"])

            newWavelengths.extend(molecule.convWavelength)
            newIntensities.extend(molecule.convIntensity)

        spectrum.append({"molecule": m.molecule, "wavelengths": newWavelengths, "intensities": newIntensities, "numberdensity": numberDensities[idx]})
    
    # fig = figure(figsize=(20,10))
    # frame = fig.add_subplot()

    # for i,m in enumerate(spectrum):
    #     frame.plot(m["wavelengths"], m["intensities"], label=f"{m['molecule']} - {numberDensities[i]:.2e} cm$^{-2}$")

    # frame.set_title(f"Slab model spectrum at temperature {temperature}K")
    # frame.set_xlabel("Wavelength [$\mu$m]")
    # frame.set_ylabel("Line intensity [mJy]")
    # frame.legend()
    # fig.tight_layout()

    # if saveFig:
    #     filename = f"slabModel_{temperature}K.pdf"
    #     savefig(filename)

    # filename = f"temperatureGrid/slabModel_{temperature}K.json"
    # saveSpectrum(filename, spectrum)

    return spectrum
    
def saveSpectrum(filename, spectrum):
    with open(filename, "w") as outfile: 
        json.dump(spectrum, outfile)

if __name__ == "__main__":
    molecule_list = ["H2O", "OH", "CO", "CO2", "HCN", "C2H2", "NH3", "SO2"]
    number_densities = [1e18, 1e15, 1e18, 1e14, 1e14, 1e14, 1e17, 1e14]
    
    # gridData("temperatureGrid", 100, 1500, 100, molecule_list, number_densities, verbose=True)

    for i in range(2,16):
        temp = i*100
        print(i)
        generateSlabModel(molecule_list, temp, number_densities, verbose=True)

    # m = Molecule("H2O", 400)
    # m.get_moldata(5,25)
    # m.generateSlab(1e18, 5, 25)

    # convolver = Convolver()
    # convolver.readResolvingModel("convolver_data/JWST_MIRI_MRS.json", wavelength_overlap="optimal")
    # # Maybe for the report, plot a mosaic of different convolver settings?
    # m = convolver.convolveData(m, 5, 25)

    # m.saveSlab()
