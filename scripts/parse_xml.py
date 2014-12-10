import pandas as pd
import glob, os
from thermopyl import Parser

XML_PATH = os.path.join(os.environ["HOME"], "dat/thermo")

data = []
compound_dict = {}
for filename in glob.glob("%s/*/*.xml" % XML_PATH):
    try:
        parser = Parser(filename)
        current_data = parser.parse()
        data.extend(current_data)
        compound_dict.update(parser.compound_name_to_formula)
    except IOError:
        continue

data = pd.DataFrame(data)
compound_dict = pd.Series(compound_dict)

data.to_hdf("./data.h5", 'data')
compound_dict.to_hdf("./compound_name_to_formula.h5", 'data')
