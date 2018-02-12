from __future__ import print_function
from __future__ import absolute_import

from .sysinfo import *
from .version import version as __version__

import os
here = os.path.abspath(os.path.dirname(__file__))
os.environ["OPENMS_DATA_PATH"] = os.path.join(here, "share/OpenMS")

try:
    from .all_modules import *
    from .python_extras import *
except Exception as e:
    print("\n")
    print("="*70)
    print("\n")
    print("maybe you miss some libraries. please run ldd (on linux) or")
    print("dependency walker (on windows) on ")
    print("\n")
    print(os.path.join(here, "pyopenms.so"))
    print("\n")
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
