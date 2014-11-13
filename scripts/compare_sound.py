import pymbar
import os
import re
from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import pandas as pd
import glob

pressure = 101.325
expt = pd.read_csv("/home/kyleb/dat/thermo/data.csv")

filenames = glob.glob("/home/kyleb/dat/TBV/*.dat")
theory = []
for filename in filenames:
    split = filename.split("_")
    cas = os.path.split(split[0])[1]
    temperature = float(split[2])
    x = pd.read_csv(filename)
    rho = x['Density (g/mL)'] * 1000.  # SI UNITS NOW
    t, g, N = pymbar.timeseries.detectEquilibration(rho)
    mu = rho[t:].mean()
    sigma = rho[t:].std()
    theory.append({"cas":cas, "temperature":temperature, "mu":mu, "sigma":sigma})

theory = pd.DataFrame(theory)
