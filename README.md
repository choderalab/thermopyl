[![Build Status](https://travis-ci.org/choderalab/thermopyl.png?branch=master)](https://travis-ci.org/choderalab/thermopyl)
[![Binstar Badge](https://binstar.org/choderalab/thermopyl-dev/badges/version.svg)](https://binstar.org/choderalab/thermopyl-dev)

ThermoPyL
=========

Tools for exploring and using the [ThermoML Archive](http://trc.nist.gov/ThermoML.html) from the [NIST TRC](http://trc.nist.gov).

## References

See the arXiv preprint:
> Towards Automated Benchmarking of Atomistic Forcefields: Neat Liquid Densities and Static Dielectric Constants from the ThermoML Data Archive
> Kyle A. Beauchamp, Julie M. Behr, Ariën S. Rustenburg, Christopher I. Bayly, Kenneth Kroenlein, John D. Chodera
> [arXiv:1506.00262](arXiv:1506.00262)

## Installation:

The easiest way to install ThermoPyL is via the [conda](http://conda.pydata.org/docs/) package manager, which comes with the [Anaconda Scientific Python Distribution](https://store.continuum.io/cshop/anaconda/):
```
conda config --add channels choderalab
conda install thermopyl
```

## Creating a local mirror of the ThermoML Archive

1.  Create a local mirror of the ThermoML Archive:
```
thermoml-update-mirror
```
By default, the archive is placed in `~/.thermoml/`.
You can use the environment variable `THERMOML_PATH` to store its location on your disk.
Re-running this command will update the local mirror with new data published in the ThermoML Archive RSS feeds.
2.  Run `thermoml-build-pandas` to create a pandas version of the database, saved as an HDF5 file in the archive directory.
3.  Use Pandas to query the experimental literature:
```python
import thermopyl
# Read ThermoML archive data into pandas dataframe
df = thermopyl.pandas_dataframe()
```
Datatypes are listed in columns (in addition to useful fields like `components`, `filename`, and `phase`):
```python
datatypes = list(df.columns)
```
Unavailable data is labeled as `NaN`. You can use this to extract useful data by querying on data that is not `NaN`.
For example, to extract dataframe rows with mass densities:
```python
# Extract rows with mass densities
densities = thermoml[np.isnan(thermoml['Mass density, kg/m3'])==False]
# Get a list of unique components
unique_components = set(x['components'])
```

## Updating a locally existing copy of the ThermoML Archive via the Python API

To update an existing locally saved archive via the Python API:
```
import thermopyl
thermopyl.update_archive()
```

Maintainers
-----------

* Kyle A Beaucamp (MSKCC)
* John D. Chodera (MSKCC)
* Arien Sebastian Rustenburg (MSKCC)
