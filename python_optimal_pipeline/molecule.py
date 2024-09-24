from molmass import Formula # type: ignore
import prodimopy.hitran as ht # type: ignore
import prodimopy.run_slab as runs # type: ignore
from scipy.constants import c,k,h
from scipy.constants import astronomical_unit as au
from scipy.constants import parsec as pc
import numpy as np

class Molecule:
    def __init__(self, moleculeName, isotopes=[1]):
        """Generates molecule data used for generating slab models
        
        Parameters:
        moleculeName    (string)    - Molecular formula of the desired molecule
        isotopes        (list[int]) - List of isotopes to be included. Default = [1]
        
        Returns:
        """
        self.molecule = moleculeName
        self.isotopologue = isotopes
        self.filepath = 'data/hitran/'+moleculeName+'.par'
        self.QTpath = 'data/QTpy/'

        f = Formula(moleculeName)
        self.mol_mass = f.mass

    def get_moldata(self, lower, upper):
        """Read the HITRAN file containing all the emission lines of the given molecule
        
        Parameters:
        lower   (float) - The lower wavelength limit in microns
        upper   (float) - The upper wavelength limit in microns"""
        self.mol_data = ht.read_hitran(self.filepath,
            self.molecule,
            self.isotopologue,
            lowerLam=lower, higherLam=upper)
        
        include = ['lambda','A','E_u','E_l','global_u','global_l','local_u','local_l','g_u','g_l']  # Include only these columns, rest are irrelevant for this project
        self.mol_data[include]

    def generateSlab(self, column_density, temperature, lower_wavelength, upper_wavelength):
        """Generates a 0D slab model using prodimopy for the given molecule at various parameters. Saves model in models folder
        
        Parameters:
        column_density      (float) - The gas column density in molecules per cm^2
        temperature         (float)     - Temperature in Kelvin of the molecule
        lower_wavelength    (float) - The lower wavelength from where to run the slab model
        upper_wavelength    (float) - The upper wavelength from where to run the slab model"""
        
        self.data = runs.run_0D_slab(Ng         = column_density,            # The gas column density in molecules per cm2
                                Tg         = temperature,               # The gas temperature
                                vturb      = 2.0,                       # Delta nu_D or the width parameter in km/s
                                molecule   = self.molecule,             # name of the molecule
                                mol_mass   = self.mol_mass,             # molecular mass in amu
                                HITRANfile = self.filepath,             # Path of the HITRAN data file
                                QTpath     = self.QTpath,               # path of the partition sums
                                isotopolog = self.isotopologue,         # list of isotopologues to include      
                                wave_mol   = [lower_wavelength,upper_wavelength],         # wavelength region in which the lines are to be picked up
                                mode       = 'line_by_line',            # "line-by-line" calculation, or include mutual "overlap"
                                output     = 'return',                  # "return" data or write to "file" or do "both"
                                )
        
        # print(dir(self.data))
        print(self.data.linedata)
        print(self.data.nlines)
        
    # def saveSlab(self):
    #     solid_angle = np.pi*(0.125*au)**2/(155*pc)**2   # Solid angle = emitting area / distance^2
    #     convWavelengths = self.convWavelength
    #     convIntensities = self.convIntensity*1e26*solid_angle
    #     with open(self.output_filename, "w") as f:
    #         f.write("convWavelength[microns],convLine_intensity[mJy]\n")
    #         for cw,ci in zip(convWavelengths, convIntensities):
    #             f.write(f"{cw},{ci}\n")


        