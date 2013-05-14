import sysinfo

from version import version as __version__

import os
here = os.path.abspath(os.path.dirname(__file__))
os.environ["OPENMS_DATA_PATH"] = os.path.join(here, "share/OpenMS")

import sys
if sys.platform != "win32":
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
    print "="*78
    print 
    print """maybe you miss some libraries. please run ldd (on linux) or dependency"""
    print "explorer (on windows) on "
    print
    print os.path.join(here, "pyopenms.so")
    print
    print "="*78
    print 
    raise e

del os, here, sys

