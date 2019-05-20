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
from __future__ import absolute_import

from .sysinfo import *
from .version import version as __version__

import os
here = os.path.abspath(os.path.dirname(__file__))
os.environ["OPENMS_DATA_PATH"] = os.path.join(here, "share/OpenMS")

import sys
if sys.platform.startswith("linux"):
    # load local shared libries before we import pyopenms.so, else
    # those are not found. setting LD_LIBRARY_PATH does not work,
    # see: http://stackoverflow.com/questions/1178094
    import ctypes
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenSwathAlgo.so"))
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenMS.so"))
    ctypes.cdll.LoadLibrary(os.path.join(here, "libSuperHirn.so"))

try:
    from .all_modules import *
    from .python_extras import *
except Exception as e:
    print("")
    print("="*70)
    print("Error when loading pyOpenMS libraries!")
    print("Libraries could not be found / could not be loaded.")
    print("")
    print("Note: when using the Spyder IDE, this error may be triggered when")
    print("the 'Automatic' backend is used. Please change this in Tools ->")
    print("Preferences -> IPython -> Graphics to 'Inline'.")
    print("")
    print("To debug this error, please run ldd (on linux) or dependency walker (on windows) on ")
    print("")
    print(os.path.join(here, "pyopenms.so"))
    print("")
    print("="*70)

    try:
        import PyQt4.QtCore
    except:
        pass
    else:
        from .qt_version_info import info

        info = "\n    ".join(info.split("\n"))

        print("""When building pyopenms qmake said:

    %s
PYQT has version %s

Maybe this causes a conflict.  You can test this by importing pyopenms
first and then import PyQt4.QtCore.
        """ % (info, PyQt4.QtCore.PYQT_VERSION_STR) )
    print("="*70)
    print("\n")
    raise e

del os, here, sys
