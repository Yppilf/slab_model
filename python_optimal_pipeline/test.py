import numpy as np
from matplotlib.pyplot import figure,savefig

loaded = np.load("spectra/spectrum0.npz")
wavelength = loaded["wavelengths"]
intensities = loaded["intensities"]
molecules = loaded["molecules"]
densities = loaded["numberDensitites"]

fig=figure()
frame = fig.add_subplot()
frame.plot(wavelength, intensities)
frame.set_xlabel("wavelength")
frame.set_ylabel("intensity")
print(molecules, densities)
savefig("spectra/spectrum0.png")