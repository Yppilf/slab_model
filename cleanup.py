import os, re

folder = "/scratch/s4950836"
newFolder = "/scratch/s4950836/spectra"
pattern = "([A-Za-z0-9]+)_([0-9]+)_([0-9.]+)_[a-z]+.npz"
pattern2 = "spectrum[0-9]+.npz"
tempSort = {}
for file in os.listdir(folder):
    match = re.match(pattern, file)
    if match:
        os.remove(f"{folder}/{file}")
        continue
    match2 = re.match(pattern2,file)
    if match2:
        os.rename(f"{folder}/{file}", f"{newFolder}/{file}")


