from convolver import Convolver
import numpy as np
from matplotlib.pyplot import figure,show,savefig,cm

convolver_settings = ["lower", "upper", "optimal", "minimal"]
resolving_power_settings = ["mean", "max", "min"]
color = iter(cm.gist_rainbow(np.linspace(0, 1, len(convolver_settings)*len(resolving_power_settings))))

fig = figure(figsize=(10,10))
frame = fig.add_subplot()

for i,setting in enumerate(convolver_settings):
    for j, setting2 in enumerate(resolving_power_settings):
        convolver = Convolver()
        convolver.readResolvingModel("./convolver_data/JWST_MIRI_MRS.json", wavelength_overlap=setting, resolving_power=setting2)

        resPowers = []
        xMin = []
        xMax = []
        for channel in convolver.data:
            resPowers.append(channel["resolving_power"])
            xMin.append(channel["wl"])
            xMax.append(channel["wu"])

        c = next(color)
        frame.hlines(resPowers, xMin, xMax, label=f"{setting} - {setting2}", color=c)

frame.set_title("Resolving power versus wavelength for JWST MIRI MRS")
frame.set_xlabel("Wavelength [$\mu$m]")
frame.set_ylabel("Resolving power")
frame.legend()
fig.tight_layout()

savefig("convolverSettings.pdf")