# input-encoding: latin-1

import distribute_setup
distribute_setup.use_setuptools()

# windows ?
import sys
iswin = sys.platform == "win32"

# import config
from env import  (OPEN_MS_SRC, OPEN_MS_BUILD_DIR, OPEN_MS_CONTRIB_BUILD_DIRS, QT_HEADERS_DIR,
                  QT_LIBRARY_DIR, QT_QTCORE_INCLUDE_DIR, MSVCR90DLL, MSVCP90DLL,
                  QT_QMAKE_VERSION_INFO, OPEN_MS_BUILD_TYPE, OPEN_MS_VERSION, LIBRARIES_EXTEND,
                  LIBRARY_DIRS_EXTEND, OPEN_MS_LIB, OPEN_SWATH_ALGO_LIB)

IS_DEBUG = OPEN_MS_BUILD_TYPE.upper() == "DEBUG"

if iswin and IS_DEBUG:
    raise Exception("building pyopenms on windows in debug mode not tested yet.")

# use autowrap to generate Cython and .cpp file for wrapping OpenMS:
import cPickle
import os
import glob
import shutil

j = os.path.join

src_pyopenms = j(OPEN_MS_SRC, "pyOpenMS")
extra_includes = glob.glob(src_pyopenms + "/extra_includes/*.h*")

for include in extra_includes:
    shutil.copy(include, "extra_includes/")


persisted_data_path = "include_dir.bin"
autowrap_include_dirs = cPickle.load(open(persisted_data_path, "rb"))

from setuptools import setup, Extension
import time

# create version information
ctime = os.stat("pyopenms").st_mtime
ts = time.gmtime(ctime)
timestamp = "%02d-%02d-%4d" % (ts.tm_mday, ts.tm_mon, ts.tm_year)

from version import version

# parse config

if OPEN_MS_CONTRIB_BUILD_DIRS.endswith(";"):
    OPEN_MS_CONTRIB_BUILD_DIRS = OPEN_MS_CONTRIB_BUILD_DIRS[:-1]

for OPEN_MS_CONTRIB_BUILD_DIR in OPEN_MS_CONTRIB_BUILD_DIRS.split(";"):
    if os.path.exists(os.path.join(OPEN_MS_CONTRIB_BUILD_DIR, "lib")):
        break


# Package data expected to be installed. On Linux the debian package
# contains share/ data and must be installed to get access to the OpenMS shared
# library.
#
if iswin:
    if OPEN_MS_BUILD_TYPE.upper() == "DEBUG":
        libraries = ["OpenMSd", "OpenSwathAlgod", "SuperHirnd", "xerces-c_3D", "QtCored4", "cblas_d"]
    else:
        libraries = ["OpenMS", "OpenSwathAlgo", "SuperHirn", "xerces-c_3", "QtCore4", "cblas"]

elif sys.platform == "linux2":
    libraries = ["OpenMS", "OpenSwathAlgo", "SuperHirn", "xerces-c", "QtCore"]
elif sys.platform == "darwin":
    libraries = ["OpenMS", "OpenSwathAlgo", "SuperHirn", "xerces-c"]
else:
    print
    print "platform ", sys.platform, "not supported yet"
    print
    exit()

library_dirs = [OPEN_MS_BUILD_DIR,
                j(OPEN_MS_BUILD_DIR, "lib"),
                j(OPEN_MS_BUILD_DIR, "bin"),
                j(OPEN_MS_BUILD_DIR, "bin", "Release"),
                j(OPEN_MS_BUILD_DIR, "Release"),
                j(OPEN_MS_CONTRIB_BUILD_DIR, "lib"),
                QT_LIBRARY_DIR,
                ]

import numpy

include_dirs = [
    "extra_includes",
    QT_HEADERS_DIR,
    QT_QTCORE_INCLUDE_DIR,
    j(OPEN_MS_CONTRIB_BUILD_DIR, "include"),
    j(OPEN_MS_CONTRIB_BUILD_DIR, "include", "libsvm"),
    # j(OPEN_MS_CONTRIB_BUILD_DIR, "src", "boost_1_52_0")
    j(OPEN_MS_CONTRIB_BUILD_DIR, "include", "boost"),
    j(OPEN_MS_CONTRIB_BUILD_DIR, "include", "WildMagic"),
    j(OPEN_MS_CONTRIB_BUILD_DIR, "include", "eigen3"),
    j(OPEN_MS_BUILD_DIR, "src/openswathalgo/include"),
    j(OPEN_MS_BUILD_DIR, "src/openms/include"),
    j(OPEN_MS_BUILD_DIR, "src/openms_gui/include"),
    j(OPEN_MS_BUILD_DIR, "src/superhirn/include"),
    j(OPEN_MS_SRC, "src/openswathalgo/include"),
    j(OPEN_MS_SRC, "src/openms/include"),
    j(OPEN_MS_SRC, "src/openms_gui/include"),
    j(OPEN_MS_SRC, "src/superhirn/include"),
    j(numpy.core.__path__[0], "include"),
]

include_dirs.extend(LIBRARIES_EXTEND)
libraries.extend(LIBRARIES_EXTEND)
library_dirs.extend(LIBRARY_DIRS_EXTEND)

extra_link_args = []

if sys.platform == "linux2":
    extra_link_args = ["-Wl,-s"]
elif sys.platform == "darwin":
    # we need to manually link to the Qt Frameworks
    extra_link_args = ["-Wl,-s",
                       "-F" + QT_LIBRARY_DIR,
                       "-framework AGL",
                       "-framework CoreFoundation",
                       "-framework CoreServices",
                       "-framework QtSvg",
                       "-framework QtWebKit",
                       "-framework QtXmlPatterns",
                       "-framework QtGui",
                       "-framework QtTest",
                       "-framework QtXml",
                       "-framework QtSql",
                       "-framework QtNetwork",
                       "-framework QtCore"]

ext = Extension(
    "pyopenms",
    sources=["pyopenms/pyopenms.cpp"],
    language="c++",
    library_dirs=library_dirs,
    libraries=libraries,
    include_dirs=include_dirs + autowrap_include_dirs,

    # /EHs is important. It sets _CPPUNWIND which causes boost to
    # set BOOST_NO_EXCEPTION in <boost/config/compiler/visualc.hpp>
    # such that  boost::throw_excption() is declared but not implemented.
    # The linker does not like that very much ...
    extra_compile_args=iswin and ["/EHs", "/bigobj"] or (IS_DEBUG and ["-g2"] or []),
    extra_link_args=extra_link_args
)


share_data = []
if iswin:
    share_data += [MSVCR90DLL, "xerces-c_3_1.dll"]

share_data.append("License.txt")

# enforce 64bit-only build as OpenMS is not available in 32bit on osx
if sys.platform == "darwin":
    os.environ['ARCHFLAGS'] = "-arch x86_64"

setup(

    name="pyopenms",
    packages=["pyopenms"],
    ext_package="pyopenms",

    version=version + "-" + timestamp,

    maintainer="Uwe Schmitt",
    maintainer_email="uschmitt@mineway.de",
    license="http://opensource.org/licenses/BSD-3-Clause",
    platforms=["any"],
    description="Python wrapper for C++ LCMS library OpenMS",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    long_description=open("README.rst").read(),
    zip_safe=False,

    url="http://open-ms.de",

    author="Uwe Schmitt",
    author_email="uschmitt@mineway.de",

    ext_modules=[ext],
    include_package_data=True  # see MANIFEST.in
)
