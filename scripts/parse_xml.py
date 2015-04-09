import pandas as pd
import glob, os
from thermopyl import Parser

XML_PATH = os.environ["THERMOML_PATH"]
filenames = glob.glob("%s/ThermoML/*/*.xml" % XML_PATH)

data = []
compound_dict = {}
for filename in filenames:
    print(filename)
    try:
        parser = Parser(filename)
        current_data = parser.parse()
        current_data = pd.DataFrame(current_data)
        data.append(current_data)
        compound_dict.update(parser.compound_name_to_formula)
    except Exception as e:
        print(e)

data = pd.concat(data, copy=False)  # Because the input data is a list of DataFrames, this saves a LOT of memory!
data.to_hdf("./data.h5", 'data')

compound_dict = pd.Series(compound_dict)
compound_dict.to_hdf("./compound_name_to_formula.h5", 'data')
