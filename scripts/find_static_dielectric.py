import re
from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import pandas as pd
import glob
import thermoml_lib

data = pd.read_hdf("./data.h5", 'data')

#experiments = ["Mass density, kg/m3"]
#experiments = ["Relative permittivity at zero frequency"]
#experiments = ["Isobaric coefficient of expansion, 1/K"]  # Essentially None
#experiments = ["Isothermal compressibility, 1/kPa"]  # Essentially None
experiments = ["Mass density, kg/m3", "Relative permittivity at zero frequency"]

ind_list = [data[exp].dropna() for exp in experiments]
ind = reduce(lambda x,y: x.index.union(y.index), ind_list)
X = data.ix[ind]

name_to_formula = pd.read_hdf("./compound_name_to_formula.h5", 'data')
name_to_formula = name_to_formula.dropna()


which_atoms = ["H","N","C","O","S","Cl","Br"]
X = X[X["Temperature, K"] > 270]
X = X[X["Temperature, K"] < 330]
X = X[X["Pressure, kPa"] > 100.]
X = X[X["Pressure, kPa"] < 102.]
X.dropna(axis=1, how='all', inplace=True)

X_is_good = {}
for k, row in X.iterrows():
    chemical_string = row.components
    chemicals = chemical_string.split("__")
    try:
        X_is_good[k] = all([thermoml_lib.is_good(name_to_formula[chemical], good_atoms=which_atoms) for chemical in chemicals])
    except KeyError:
        print("Warning, could not find %d %s" % (k, chemical_string))
        X_is_good[k] = False

X_is_good = pd.Series(X_is_good)
X["is_good"] = X_is_good
X = X[X.is_good]

X["n_components"] = X.components.apply(lambda x: len(x.split("__")))
X = X[X.n_components == 1]
X.dropna(axis=1, how='all', inplace=True)

X["n_heavy_atoms"] = X.components.apply(lambda x: thermoml_lib.count_atoms(name_to_formula[x]))
X = X[X.n_heavy_atoms <= 10]
X.dropna(axis=1, how='all', inplace=True)

X["n_atoms"] = X.components.apply(lambda x: thermoml_lib.count_atoms(name_to_formula[x], which_atoms=which_atoms))
X = X[X.n_atoms <= 100]
X.dropna(axis=1, how='all', inplace=True)

X["smiles"] = X.components.apply(lambda x: cirpy.resolve(x, "smiles"))  # This should be cached via sklearn.
X = X[X.smiles != None]
X = X.ix[X.smiles.dropna().index]

    
X["cas"] = X.components.apply(lambda x: thermoml_lib.get_first_entry(cirpy.resolve(x, "cas")))  # This should be cached via sklearn.
X = X[X.cas != None]
X = X.ix[X.cas.dropna().index]

X["Pressure, kPa"] = 101.325  # Assume everything within range is comparable.  

# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!
bad_filenames = ["./10.1016/j.fluid.2013.12.014.xml"]
X = X[~X.filename.isin(bad_filenames)]
# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!# SEE GOOGLE DOC!!!!!!

mu = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].mean()
sigma = X.groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa"])[experiments].std().dropna()

mu = mu.dropna()

q = mu.reset_index()
len(q.components.unique())

# 400 measurements, 48 components.

q.to_csv("./data_dielectric.csv")


#Y = X[X.components == u'2-propanol'].groupby(["components", "smiles", "cas", "Temperature, K", "Pressure, kPa", "filename"])[experiments].mean().dropna()
#Y = X[X.components == u'2-propanol'][experiments[1:]].dropna()
