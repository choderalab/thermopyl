#!/usr/bin/env python
"""
Parse ThermoML XML files in the local ThermoML Archive mirror.

"""
import pandas as pd
import glob, os, os.path
from thermopyl import Parser

def main():
    try:
        XML_PATH = os.environ["THERMOML_PATH"]
    except:
        XML_PATH = os.path.join(os.environ["HOME"], '.thermoml')

    filenames = glob.glob("%s/*.xml" % XML_PATH)

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

    data = pd.concat(data, copy=False, ignore_index=True)  # Because the input data is a list of DataFrames, this saves a LOT of memory!  Ignore the index to return unique index.
    data.to_hdf("%s/data.h5" % XML_PATH, 'data')

    compound_dict = pd.Series(compound_dict)
    compound_dict.to_hdf("%s/compound_name_to_formula.h5" % XML_PATH, 'data')

    return

if __name__ == '__main__':
    main()
