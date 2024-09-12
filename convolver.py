import json
from scipy.constants import astronomical_unit as au
from scipy.constants import parsec as pc
import numpy as np

class Convolver:
    def __init__(self):
        pass

    def readResolvingModel(self, filename, wavelength_overlap="lower", resolving_power="mean"):
        """Read a JSON resolving model into memory
        
        Parameters:
        filename            (string)    - The filepath to the JSON file containing the data
        wavelength_overlap  (string)    - What method to apply for overlapping wavelength ranges of channels. Default = lower
            Options:
                lower   - Always takes the maximum of the lower channel range, continuing the second channel from that point
                upper   - Always takes the maximum of the higher channel range, using the lower channel until that point
                optimal - Takes the wavelength range of the channel with highest resolving power, cutting the other channel at that point
                minimal - Takes the wavelength range of the channel with the lowest resolving power, cutting the other channel at that point
        resolving_power     (string)    - The method of determining the resolving power used. Default = mean
            Options:
                mean    - Takes the average between upper and lower limit of resolving power of the channel
                max     - Takes the upper limit
                min     - Takes the lower limit"""


        # To whomever might need to read this, sorry
        with open(filename) as json_data:
            d = json.load(json_data)

        # Determine used resolving power
        for channel in d:
            u = channel["resolving_power_upper"]
            l = channel["resolving_power_lower"]

            if resolving_power=="mean":
                channel["resolving_power"] = (u+l)/2
            elif resolving_power=="max":
                channel["resolving_power"] = u
            elif resolving_power=="min":
                channel["resolving_power"] = l

        # Determine used wavelength range
        for i,channel in enumerate(d):
            if wavelength_overlap=="lower":
                if i==0:
                    channel["wl"] = channel["wavelength_lower"]
                else:
                    channel["wl"] = d[i-1]["wu"]
                channel["wu"] = channel["wavelength_upper"]
            elif wavelength_overlap=="upper":
                if i==len(d)-1:
                    channel["wu"] = channel["wavelength_upper"]
                else:
                    channel["wu"] = d[i+1]["wavelength_lower"]
                channel["wl"] = channel["wavelength_lower"]
            elif wavelength_overlap=="optimal":
                if i==0:
                    channel["wl"] = channel["wavelength_lower"]
                    if channel["resolving_power"] > d[i+1]["resolving_power"]:
                        channel["wu"] = channel["wavelength_upper"]
                    else:
                        channel["wu"] = d[i+1]["wavelength_lower"]
                elif i==len(d)-1:
                    channel["wu"] = channel["wavelength_upper"]
                    if channel["resolving_power"] > d[i-1]["resolving_power"]:
                        channel["wl"] = channel["wavelength_lower"]
                    else:
                        channel["wl"] = d[i-1]["wavelength_upper"]
                else:
                    if channel["resolving_power"] > d[i-1]["resolving_power"]:
                        channel["wl"] = channel["wavelength_lower"]
                    else:
                        channel["wl"] = d[i-1]["wavelength_upper"]

                    if channel["resolving_power"] > d[i+1]["resolving_power"]:
                        channel["wu"] = channel["wavelength_upper"]
                    else:
                        channel["wu"] = d[i+1]["wavelength_lower"]
            elif wavelength_overlap=="minimal":
                if i==0:
                    channel["wl"] = channel["wavelength_lower"]
                    if channel["resolving_power"] < d[i+1]["resolving_power"]:
                        channel["wu"] = channel["wavelength_upper"]
                    else:
                        channel["wu"] = d[i+1]["wavelength_lower"]
                elif i==len(d)-1:
                    channel["wu"] = channel["wavelength_upper"]
                    if channel["resolving_power"] < d[i-1]["resolving_power"]:
                        channel["wl"] = channel["wavelength_lower"]
                    else:
                        channel["wl"] = d[i-1]["wavelength_upper"]
                else:
                    if channel["resolving_power"] < d[i-1]["resolving_power"]:
                        channel["wl"] = channel["wavelength_lower"]
                    else:
                        channel["wl"] = d[i-1]["wavelength_upper"]

                    if channel["resolving_power"] < d[i+1]["resolving_power"]:
                        channel["wu"] = channel["wavelength_upper"]
                    else:
                        channel["wu"] = d[i+1]["wavelength_lower"]

        self.data = d
                
    def convolveData(self, molecule, lower, upper):
        """Applies resolving power calculations on the generated spectrum by the molecule. Assumes the entire range of the molecule to be within one channel
        
        Parameters:
        molecule    (Molecule)  - Molecule object for which to convolve the data
        lower       (float)     - Lower wavelength for convolving in microns
        upper       (float)     - Upper wavelength for convolving in microns"""

        data = molecule.data

        for channel in self.data:
            if lower >= channel["wl"] and lower < channel["wu"]:
                res = channel["resolving_power"]
                break

        data.convolve(R=res,lambda_0=lower,lambda_n=upper,verbose=False)
        molecule.convWavelength = data.convWavelength
        solid_angle = np.pi*(0.125*au)**2/(155*pc)**2   # Solid angle = emitting area / distance^2
        molecule.convIntensity = data.convLTEflux*1e26*solid_angle
        return molecule