[![Build Status](https://travis-ci.org/choderalab/thermopyl.png)](https://travis-ci.org/choderalab/thermopyl)
[![Binstar Badge](https://binstar.org/omnia/thermopyl/badges/version.svg)](https://binstar.org/omnia/thermopyl)

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
conda config --add channels http://conda.binstar.org/omnia
conda install thermopyl
```

## Creating a local mirror of the ThermoML Archive
1.  Obtain an archive of the ThermoML archive:
```
thermoml-archive-update
```
By default, the archive is placed in `~/.thermoml/`.
You can use the environment variable `THERMOML_PATH` to store its location on your disk.

2.  Run `thermoml-to-pandas` to create a pandas version of the database, saved as an HDF5 file

3.  Use Pandas to query the experimental literature

## Updating a locally existing copy of the ThermoML Archive

To update an existing locally saved archive:
```
import thermopyl
thermopyl.update_archive()
```

Maintainers
-----------

* John D. Chodera (MSKCC)
* Patrick B. Grinaway (MSKCC)
* Arien Sebastian Rustenburg (MSKCC)

Contributors
------------

* Kyle A Beaucamp (MSKCC)
