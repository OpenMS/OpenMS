#!/usr/bin/python
# -*- encoding: utf8 -*-
"""Python bindings to the OpenMS C++ library.

The pyOpenMS package contains Python bindings for a large part of the OpenMS
library (http://www.open-ms.de) for mass spectrometry based proteomics. It thus
provides providing facile access to a feature-rich, open-source algorithm
library for mass-spectrometry based proteomics analysis. These Python bindings
allow raw access to the data-structures and algorithms implemented in OpenMS,
specifically those for file access (mzXML, mzML, TraML, mzIdentML among
others), basic signal processing (smoothing, filtering, de-isotoping and
peak-picking) and complex data analysis (including label-free, SILAC, iTRAQ and
SWATH analysis tools).

For further documentation, please see https://pyopenms.readthedocs.io

Please cite:

    Röst HL, Schmitt U, Aebersold R, Malmström L.
    pyOpenMS: a Python-based interface to the OpenMS mass-spectrometry algorithm library.
    Proteomics. 2014 Jan;14(1):74-7. doi: 10.1002/pmic.201300246.

"""
from __future__ import print_function

from ._sysinfo import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
from ._version import version as __version__

import os
here = os.path.abspath(os.path.dirname(__file__))

default_openms_data_path = os.path.join(here, "share/OpenMS")
env_openms_data_path = os.environ.get("OPENMS_DATA_PATH")

if os.path.exists(default_openms_data_path):
    if not env_openms_data_path:
        os.environ["OPENMS_DATA_PATH"] = default_openms_data_path
    else:
        print(
            "Warning: OPENMS_DATA_PATH environment variable already exists. "
            "pyOpenMS will use it ({env}) to locate data in the OpenMS share folder "
            "(e.g., the unimod database), instead of the default ({default})."
            .format(env=env_openms_data_path, default=default_openms_data_path)
        )
else:
    if not env_openms_data_path:
         print(
            "Warning: OPENMS_DATA_PATH environment variable not found and no share directory was installed. "
            "Some functionality might not work as expected."
        )

import sys
# on conda the libs will be installed to the general conda lib path which is available during load.
# try to skip this loading if we do not ship the libraries in the package (e.g. as wheel via pip)
# TODO check if this can be completely removed by now or e.g. by baking in an RPATH into the pyopenms*.so's
if sys.platform.startswith("linux") and os.path.exists(os.path.join(here, "libOpenMS.so")):
    # load local shared libraries before we import pyopenms*.so, else
    # those are not found. setting LD_LIBRARY_PATH does not work,
    # see: http://stackoverflow.com/questions/1178094
    import ctypes
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenSwathAlgo.so"))
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenMS.so"))

try:
    from ._all_modules import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
    from ._python_extras import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
    # This has to be imported after all_modules so it can augment the core datastructures with dataframe
    # export capabilities
    from ._dataframes import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
except Exception as e:
    print("")
    print("="*70)
    print("Error when loading pyOpenMS libraries!")
    print("Libraries could not be found / could not be loaded.")
    print("")
    print("To debug this error, please run ldd (on linux), otool -L (on macOS) or dependency walker (on windows) on ")
    print("")
    print(os.path.join(here, "pyopenms*.so"))
    print("")
    print("="*70)

    try:
        import PyQt5.QtCore
    except:
        pass
    else:
        from ._qt_version_info import info

        info = "\n    ".join(info.split("\n"))

        print("""PyQt5 was found to be installed. When building pyopenms qmake said:

    %s
PYQT has version %s

This might cause a conflict if both are loaded.  You can test this by importing pyopenms
first and then import PyQt5.QtCore.
        """ % (info, PyQt5.QtCore.PYQT_VERSION_STR) )
        
        print("Note: when using the Spyder IDE, the usage of PyQt might be circumvented")
        print("by not using the 'Automatic' backend. Please change this in Tools ->")
        print("Preferences -> IPython -> Graphics to 'Inline'.")
        print("In general, try to install everything with conda in the same environment to make sure Qt is used in the same version.")
        print("")
        print("="*70)
        print("\n")
    raise e

del os, here, sys
