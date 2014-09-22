import pandas as pd
import glob
import thermoml_lib

data = []
compound_dict = {}
for filename in glob.glob("./*/*.xml")[0:250]:
    try:
        parser = thermoml_lib.Parser(filename)
        current_data = parser.parse()
        data.extend(current_data)
        compound_dict.update(parser.compound_name_to_formula)
    except IOError:
        continue

data = pd.DataFrame(data)
data.to_hdf("./data.h5", 'data')

compound_dict = pd.Series(compound_dict)
compound_dict.to_hdf("./compound_name_to_formula.h5", 'data')

#chemicals = set()
#for chemlist in data.components.apply(lambda x: x.split("__")).values:
#    for chem in chemlist:
#        chemicals.add(chem)
