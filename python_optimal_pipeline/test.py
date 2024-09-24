import numpy as np
from matplotlib.pyplot import figure,savefig

loaded = np.load("H2O_1300_18.0_lower.npz")
wavelength = loaded["wavelengths"]
intensities = loaded["intensities"]
# molecules = loaded["molecules"]
# densities = loaded["numberDensitites"]
print(wavelength)
print(max(intensities))
# print(molecules)
# print(densities)
# fig=figure()
# frame = fig.add_subplot()
# frame.plot(wavelength, intensities)
# savefig("spectrum0.png")