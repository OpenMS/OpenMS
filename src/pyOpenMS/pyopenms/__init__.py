#!/usr/bin/python
# -*- encoding: utf8 -*-
"""Python bindings to the OpenMS C++ library.

The pyOpenMS package contains Python bindings for a large part of the OpenMS
library (https://openms.de) for mass spectrometry based proteomics. It thus
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
import os
import sys
import warnings

from ._sysinfo import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)

try:
    import importlib.metadata
    __version__ = importlib.metadata.version("pyopenms")
except Exception:
    __version__ = "0+unknown"

here = os.path.abspath(os.path.dirname(__file__))

default_openms_data_path = os.path.join(here, "share", "OpenMS")
env_openms_data_path = os.environ.get("OPENMS_DATA_PATH")

if os.path.exists(default_openms_data_path):
    if not env_openms_data_path:
        os.environ["OPENMS_DATA_PATH"] = default_openms_data_path
    elif os.path.abspath(env_openms_data_path) != os.path.abspath(default_openms_data_path):
        warnings.warn(
            "Warning: OPENMS_DATA_PATH environment exists and points to a different location "
            "than the default share directory. "
            "pyOpenMS will use it ({env}) to locate data in the OpenMS share folder "
            "(e.g., the unimod database), instead of the default ({default})."
            .format(env=env_openms_data_path, default=default_openms_data_path)
        )
else:
    if not env_openms_data_path:
        warnings.warn(
            "Warning: OPENMS_DATA_PATH environment variable not found and no share directory was installed. "
            "Some functionality might not work as expected."
        )

# Pre-register pyarrow filesystem handlers before C++ Arrow loads.
# This prevents "Attempted to register factory for scheme 'file'" errors
# when both pyarrow (Python) and Arrow C++ (via OpenMS WITH_PARQUET) are present.
# See: https://github.com/apache/arrow/issues/44696
try:
    from pyarrow import fs as _pa_fs
    _pa_fs.LocalFileSystem()  # Force filesystem handler registration
    del _pa_fs
except ImportError:
    pass  # pyarrow not installed, no conflict possible

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
    # This has to be imported after all_modules so it can augment the core datastructures with dataframe
    # export capabilities
    from ._dataframes import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
    from ._python_extras import *  # pylint: disable=wildcard-import; lgtm(py/polluting-import)
except Exception as e:
    print(f"""
======================================================================
Error when loading pyOpenMS libraries!
Libraries could not be found / could not be loaded.

To debug this error, please run ldd (on linux), otool -L (on macOS) or dependency walker (on windows) on

{os.path.join(here, "pyopenms*.so")}

======================================================================
""")

    try:
        import PyQt6.QtCore
    except:
        pass
    else:
        from ._dependency_version_info import qt_version

        info = "\n    ".join(qt_version.split("\n"))

        warnings.warn(
            f"""PyQt6 was found to be installed.
    pyopenms was built with Qt version: {info}
    PyQt6 version detected: {PyQt6.QtCore.PYQT_VERSION_STR}

    This may cause a conflict if both are loaded. To test for issues, try importing pyopenms
    first, then import PyQt6.QtCore.

    Note: If you are using the Spyder IDE, you can avoid PyQt conflicts by setting
    the graphics backend to 'Inline' (Tools → Preferences → IPython Console → Graphics).
    In general, ensure all dependencies are installed within the same environment (e.g., via conda)
    to guarantee compatible Qt versions.

    {"="*70}
    """
        )
    raise e

del os, here, sys
