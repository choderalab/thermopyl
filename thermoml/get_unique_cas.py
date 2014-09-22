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
