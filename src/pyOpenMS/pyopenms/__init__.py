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

from .sysinfo import *
from .version import version as __version__

import os
here = os.path.abspath(os.path.dirname(__file__))

default_openms_data_path = os.path.join(here, "share/OpenMS")
env_openms_data_path = os.environ.get("OPENMS_DATA_PATH")

if not env_openms_data_path:
    os.environ["OPENMS_DATA_PATH"] = default_openms_data_path
else:
    print(
        "Warning: OPENMS_DATA_PATH environment variable already exists. "
        "pyOpenMS will use it ({env}) to locate data in the OpenMS share folder "
        "(e.g., the unimod database), instead of the default ({default})."
        .format(env=env_openms_data_path, default=default_openms_data_path)
    )

try:
    from .all_modules import *
    from .python_extras import *
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
        from .qt_version_info import info

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
