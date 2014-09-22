import time
import re
from rdkit import Chem
from rdkit.Chem import AllChem
import cirpy
import pandas as pd
import glob
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

u = X.components.unique()
mapping = {}
for key in u:
    smiles = cirpy.resolve(key, "smiles")
    mapping[key] = smiles
    print(key, smiles)
    time.sleep(1)

mapping = pd.Series(mapping)
mapping.to_hdf("./name_mapping.h5", 'data')

mapping = pd.read_hdf("./name_mapping.h5", "data")

to_cas = {}
for key in mapping.index:
    cas = cirpy.resolve(key, "cas")
    if type(cas) == type([]):
        cas = cas[0]
    to_cas[key] = cas
    print(key, cas)
    time.sleep(1)

to_cas = pd.Series(to_cas)
to_cas.to_hdf("./cas_mapping.h5", 'data')
