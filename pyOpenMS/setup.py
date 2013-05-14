#input-encoding: latin-1

import distribute_setup
distribute_setup.use_setuptools()

# windows ?
import sys
iswin = sys.platform == "win32"

# import config
from env import *
IS_DEBUG =  OPEN_MS_BUILD_TYPE.upper() == "DEBUG"

if iswin and IS_DEBUG:
    raise Exception("building pyopenms on windows in debug mode not tested yet.")


# use autowrap to generate cython  and cpp file for wrapping openms:
import autowrap.Main
import glob
import cPickle
import os.path

pxd_files = glob.glob("pxds/*.pxd")
addons = glob.glob("addons/*.pyx")
converters = glob.glob("converters/*.py")

mtimes = [os.path.getmtime(f) for f in pxd_files + addons + converters]

if os.path.exists("pyopenms/pyopenms.pyx"):
    mtime_result = os.path.getmtime("pyopenms/pyopenms.pyx")
else:
    mtime_result = 0

persisted_data_path = "include_dir.bin"

if not os.path.exists(persisted_data_path)\
    or any(m > mtime_result for m in mtimes):

    extra_cimports = [ # "from libc.stdint cimport *",
                "from libc.stddef cimport *",
                "from UniqueIdInterface cimport setUniqueId as _setUniqueId",
                "from Map cimport Map as _Map",
                "cimport numpy as np"]
    autowrap_include_dirs = autowrap.Main.run(pxd_files,
                                            addons,
                                            converters,
                                            "pyopenms/pyopenms.pyx",
                                            extra_cimports)
    cPickle.dump(autowrap_include_dirs, open(persisted_data_path, "wb"))

else:
    autowrap_include_dirs = cPickle.load(open(persisted_data_path, "rb"))

from setuptools import setup, Extension
import os, shutil
import time

# create version information

ctime = os.stat("pyopenms").st_mtime
ts = time.gmtime(ctime)
timestamp = "%02d-%02d-%4d" % (ts.tm_mday, ts.tm_mon, ts.tm_year)

from version import version
full_version= "%s_%s" % (version, timestamp)

print >> open("pyopenms/version.py", "w"), "version=%r\n" % version


# parse config

if OPEN_MS_CONTRIB_BUILD_DIRS.endswith(";"):
    OPEN_MS_CONTRIB_BUILD_DIRS = OPEN_MS_CONTRIB_BUILD_DIRS[:-1]

for OPEN_MS_CONTRIB_BUILD_DIR in  OPEN_MS_CONTRIB_BUILD_DIRS.split(";"):
    if os.path.exists(OPEN_MS_CONTRIB_BUILD_DIR):
        break


j=os.path.join

# package data expected to be installed. on linux the debian package
# contains share/ data and must be intalled to get access to the openms shared
# library.
#
if iswin:
    shutil.copy(OPEN_MS_LIB, "pyopenms")
    shutil.copy(OPEN_SWATH_ALGO_LIB, "pyopenms")

    shutil.copy(MSVCR90DLL, "pyopenms")
    shutil.copy(MSVCP90DLL, "pyopenms")

    shutil.copy(j(OPEN_MS_CONTRIB_BUILD_DIR, "lib", "xerces-c_3_0.dll"),\
                    "pyopenms")

    if OPEN_MS_BUILD_TYPE.upper() == "DEBUG":
        libraries=["OpenMSd", "xerces-c_3D", "QtCored4", "gsl_d", "cblas_d"]
    else:
        libraries=["OpenMS", "xerces-c_3", "QtCore4", "gsl", "cblas"]

elif sys.platform == "linux2":

    libraries=["OpenMS", "OpenSwathAlgo", "xerces-c", "QtCore", "gsl",
                        "gslcblas",
              ]

    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenMS.so" ), "pyopenms")
    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenSwathAlgo.so" ), "pyopenms")

else:
    print
    print "platform ", sys.platform, "not supported yet"
    print
    exit()

library_dirs=[ OPEN_MS_BUILD_DIR,
               j(OPEN_MS_BUILD_DIR,"lib"),
               j(OPEN_MS_CONTRIB_BUILD_DIR,"lib"),
               QT_LIBRARY_DIR,
              ]


import numpy

include_dirs=[
    QT_HEADERS_DIR,
    QT_QTCORE_INCLUDE_DIR,
    j(OPEN_MS_CONTRIB_BUILD_DIR, "include"),
    j(OPEN_MS_CONTRIB_BUILD_DIR, "src", "boost_1_42_0", "include", "boost-1_42"),
    j(OPEN_MS_BUILD_DIR ,  "include"),
    j(OPEN_MS_SRC ,  "include"),
    j(numpy.core.__path__[0],"include"),
             ]



ext = Extension(
        "pyopenms",
        sources = ["pyopenms/pyopenms.cpp"],
        language="c++",
        library_dirs = library_dirs,
        libraries = libraries,
        include_dirs = include_dirs + autowrap_include_dirs,

        # /EHs is important. It sets _CPPUNWIND which causes boost to
        # set BOOST_NO_EXCEPTION in <boost/config/compiler/visualc.hpp>
        # such that  boost::throw_excption() is declared but not implemented.
        # The linker does not like that very much ...
        extra_compile_args = iswin and [ "/EHs"] or (IS_DEBUG and ["-g2"] or [])

    )


source_share = j(OPEN_MS_SRC, "share")
target_share = "pyopenms/share"
if os.path.exists(target_share):
    shutil.rmtree(target_share)
shutil.copytree(source_share, target_share, ignore=lambda *a: ["examples"])

share_data = []

if iswin:
    share_data += [MSVCR90DLL, "xerces-c_3_0.dll"]


share_data.append("License.txt")

setup(

    name = "pyopenms",
    packages = ["pyopenms"],
    ext_package = "pyopenms",

    version = full_version,

    url="http://github.com/uweschmitt/pyopenms",

    author="Uwe Schmitt",
    author_email="uschmitt@mineway.de",

    ext_modules = [ext ],

    # setup_requires=["autowrap", "cython"],

    # package_data= { "pyopenms": share_data },
    include_package_data = True # see MANIFEST.in
)
