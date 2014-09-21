import copy
import pandas as pd
import glob
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`

def parse(filename, compound_dict):
    print(filename)
    root = thermoml_schema.CreateFromDocument(open(filename).read())

    for Compound in root.Compound:
        sCommonName = Compound.sCommonName[0]
        sFormulaMolec = Compound.sFormulaMolec
        compound_dict[sCommonName] = sFormulaMolec


compound_dict = {}
for filename in glob.glob("./*/*.xml"):
    try:
        parse(filename, compound_dict)
    except IOError:
        continue

X = pd.Series(compound_dict)

X.to_hdf("./compounds.h5", 'data')
