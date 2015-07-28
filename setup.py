"""ThermoPyl: Python tools for ThermoML
"""

from __future__ import print_function, absolute_import

DOCLINES = __doc__.split("\n")

import os
import sys
import glob
import traceback
import numpy as np
from os.path import join as pjoin
from setuptools import setup, Extension, find_packages
try:
    sys.dont_write_bytecode = True
    sys.path.insert(0, '.')
    from basesetup import write_version_py, CompilerDetection, check_dependencies
finally:
    sys.dont_write_bytecode = False


if '--debug' in sys.argv:
    sys.argv.remove('--debug')
    DEBUG = True
else:
    DEBUG = False

# #########################
VERSION = '0.9'
ISRELEASED = False
__version__ = VERSION
# #########################

CLASSIFIERS = """\
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)
Programming Language :: C++
Programming Language :: Python
Development Status :: 4 - Beta
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS
Programming Language :: Python :: 2
Programming Language :: Python :: 2.6
Programming Language :: Python :: 2.7
Programming Language :: Python :: 3
Programming Language :: Python :: 3.3
Programming Language :: Python :: 3.4
"""

extensions = []

setup(name='thermopyl',
      author='Kyle Beauchamp',
      author_email='kyle.beauchamp@choderalab.org',
      description=DOCLINES[0],
      long_description="\n".join(DOCLINES[2:]),
      version=__version__,
      url='https://github.com/choderalab/thermopyl',
      platforms=['Linux', 'Mac OS-X', 'Unix'],
      classifiers=CLASSIFIERS.splitlines(),
      packages=['thermopyl', 'thermopyl.tests', 'thermopyl.scripts'],
      package_data={'thermopyl': ['data/*']},  # Install all data directories of the form /data/
      zip_safe=False,
      ext_modules=extensions,
      install_requires=[
        'six',
        'pandas',
        'pyxb==1.2.4',
        'feedparser',
        'tables',
        ],
      entry_points={'console_scripts': [
            'thermoml-update-mirror = thermopyl.scripts.update_archive:main',
            'thermoml-build-pandas = thermopyl.scripts.parse_xml:main',
            ]})
