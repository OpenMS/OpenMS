import sysinfo

from version import version as __version__

import os
here = os.path.abspath(os.path.dirname(__file__))
os.environ["OPENMS_DATA_PATH"] = os.path.join(here, "share/OpenMS")

import sys
if sys.platform == "linux2":
    # load local shared libries before we import pyopenms.so, else
    # those are not found. setting LD_LIBRARY_PATH does not work,
    # see: http://stackoverflow.com/questions/1178094
    import ctypes
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenSwathAlgo.so"))
    ctypes.cdll.LoadLibrary(os.path.join(here, "libOpenMS.so"))

try:
    from pyopenms import *
except Exception, e:
    print
    print "="*70
    print
    print "maybe you miss some libraries. please run ldd (on linux) or"
    print "dependency walker (on windows) on "
    print
    print os.path.join(here, "pyopenms.so")
    print
    try:
        import PyQt4.QtCore
    except:
        pass
    else:
        from qt_version_info import info

        info = "\n    ".join(info.split("\n"))

        print """When building pyopenms qmake said:

    %s
PYQT has version %s

Maybe this causes a conflict.  You can test this by importing pyopenms
first and then import PyQt4.QtCore.
        """ % (info, PyQt4.QtCore.PYQT_VERSION_STR)
    print "="*70
    print
    raise e

del os, here, sys

