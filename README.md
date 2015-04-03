ThermoML
========

Tools for ThermoML parsing.

Requirements:
Pandas
pyxb version 1.2.4


To use:

1.  Install the thermopyl library:
python setup.py install

2.  Obtain an archive of the ThermoML archive.  Use environment variable THERMOML_PATH to store its location on your disk.
2.  Execute parse_xml.py to create a pandas version of the database, saved as an HDF5 file
3.  Use Pandas to query the experimental literature


Usage:

To update an exting locally saved archive:

import thermopyl
thermopyl.update_archive()
