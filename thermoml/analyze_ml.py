import re
from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import pandas as pd
import glob
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

def is_good(formula_string):
    good = set(["H", "N", "C", "O", "S", ""])
    elements = set(re.findall("[A-Z][a-z]?", formula_string))    
    if good.intersection(elements) is None:  # Has no good elements
        return False
    if len(elements.difference(good)) > 0:  # Has an unwanted element
        return False
    return True

def count_heavy_atoms(formula_string):
    heavy_atoms = ["N", "C", "O", "S"]
    elements = re.findall("[A-Z][a-z]?\d?\d?\d?", formula_string)
    print(elements)
    n_heavy = 0
    for s in elements:
        try:
            n_atoms = int(re.split("[A-Z][a-z]?", s)[1])
        except:
            n_atoms = 1
        atom = re.split("\d?\d?\d?", s)[0]
        print(s, atom, n_atoms)
        if atom in heavy_atoms:
            n_heavy += n_atoms
    return n_heavy
        

data = pd.read_hdf("./data.h5", 'data')
name_to_formula = pd.read_hdf("./compounds.h5", 'data')
name_to_formula = name_to_formula.dropna()

X = data.ix[data["Mass density, kg/m3"].dropna().index]
name_to_smiles = pd.read_hdf("./name_mapping.h5", 'data')
X["smiles"] = X.components.apply(lambda x: cirpy.resolve(x, "smiles"))  # This should be cached via sklearn.
X = X[X.smiles != None]
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


mu = X.groupby(["smiles", "Temperature, K", "Pressure, kPa"])["Mass density, kg/m3"].mean()
sigma = X.groupby(["smiles", "Temperature, K", "Pressure, kPa"])["Mass density, kg/m3"].std().dropna()

