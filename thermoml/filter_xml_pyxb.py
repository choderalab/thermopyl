import copy
import pandas as pd
import glob
import thermoml_schema  # Obtained by `wget http://media.iupac.org/namespaces/ThermoML/ThermoML.xsd` and `pyxbgen ThermoML.xsd`
import thermoml_lib

data = []
for filename in glob.glob("./*/*.xml")[0:250]:
    try:
        alldata, root = parse(filename)
    except IOError:
        continue
    for d in alldata:
        if True or u'Mass density, kg/m3' in d:
            data.append(d)

data = pd.DataFrame(data)
#data.to_hdf("./data.h5", 'data')

X = data.ix[data["Mass density, kg/m3"].dropna().index]
X = X[X["Temperature, K"] > 270]
X = X[X["Temperature, K"] < 330]
X = X[X["Pressure, kPa"] > 50.]
X = X[X["Pressure, kPa"] < 150.]
X.dropna(axis=1, how='all', inplace=True)


chemicals = X.components.apply(lambda x: x.split("__")).values
s = set()
for chemlist in chemicals:
    for chem in chemlist:
        s.add(chem)
chemicals = s
