import re
from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import pandas as pd
import glob
import thermo_lib

data = pd.read_hdf("./data.h5", 'data')

experiment = "Mass density, kg/m3"

X = data.ix[data[experiment].dropna().index]

name_to_formula = pd.read_hdf("./compounds.h5", 'data')
name_to_formula = name_to_formula.dropna()

name_to_smiles = pd.read_hdf("./name_mapping.h5", 'data')

name_to_cas = pd.read_hdf("./cas_mapping.h5", 'data')
name_to_cas = name_to_cas.dropna()

X = X[X["Temperature, K"] > 270]
X = X[X["Temperature, K"] < 330]
X = X[X["Pressure, kPa"] > 50.]
X = X[X["Pressure, kPa"] < 150.]
X.dropna(axis=1, how='all', inplace=True)

X_is_good = {}
for k, row in X.iterrows():
    chemical_string = row.components
    chemicals = chemical_string.split("__")
    try:
        X_is_good[k] = all([is_good(name_to_formula[chemical]) for chemical in chemicals])
    except KeyError:
        print("Warning, could not find %d %s" % (k, chemical_string))
        X_is_good[k] = False

X_is_good = pd.Series(X_is_good)
X["is_good"] = X_is_good
X = X[X.is_good]

X["n_components"] = X.components.apply(lambda x: len(x.split("__")))
X = X[X.n_components == 1]
X.dropna(axis=1, how='all', inplace=True)

X["n_heavy_atoms"] = X.components.apply(lambda x: count_heavy_atoms(name_to_formula[x]))
X = X[X.n_heavy_atoms <= 10]
X.dropna(axis=1, how='all', inplace=True)

X["smiles"] = X.components.apply(lambda x: cirpy.resolve(x, "smiles"))  # This should be cached via sklearn.
X = X[X.smiles != None]
X = X.ix[X.smiles.dropna().index]

    
X["cas"] = X.components.apply(lambda x: first_entry(cirpy.resolve(x, "cas")))  # This should be cached via sklearn.
X = X[X.cas != None]
X = X.ix[X.cas.dropna().index]



mu = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])["Mass density, kg/m3"].mean()
sigma = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])["Mass density, kg/m3"].std().dropna()

q = mu.reset_index()
q.to_csv("./densities.csv")


