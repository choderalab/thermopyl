import pandas as pd
from simtk import unit as u
from trustbutverify import mixture_system

#data = pd.read_csv("./data_100x.csv")
data = pd.read_csv("./data_nitrobenzene.csv")

for k0, k1, components, smiles, cas, temperature, pressure, density in data.itertuples():
    print(k0, k1, components, smiles, cas, temperature, pressure, density)
    model = mixture_system.MixtureSystem([cas], [1000], temperature * u.kelvin, pressure * u.kilopascal)
    model.build()
    model.equilibrate()
    model.production()
    break
