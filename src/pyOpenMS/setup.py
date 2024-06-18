# input-encoding: latin-1
from __future__ import print_function

# windows ?
import sys
iswin = sys.platform == "win32"

# osx ?
isosx = sys.platform == "darwin"
if isosx:
    import platform
    osx_ver = platform.mac_ver()[0] #e.g. ('10.15.1', ('', '', ''), 'x86_64')

import sys
no_optimization = False

if "--no-optimization" in sys.argv:
    no_optimization = True
    sys.argv.remove("--no-optimization")

# import config
from env import  (OPEN_MS_COMPILER, PYOPENMS_SRC_DIR, OPEN_MS_GIT_BRANCH, OPEN_MS_BUILD_DIR, OPEN_MS_CONTRIB_BUILD_DIRS,
                  QT_INSTALL_LIBS, QT_INSTALL_BINS, MSVS_RTLIBS,
                  OPEN_MS_BUILD_TYPE, OPEN_MS_VERSION, INCLUDE_DIRS_EXTEND, LIBRARIES_EXTEND,
                  LIBRARY_DIRS_EXTEND, OPEN_MS_LIB, OPEN_SWATH_ALGO_LIB, PYOPENMS_INCLUDE_DIRS,
                  PY_NUM_MODULES, PY_NUM_THREADS, SYSROOT_OSX_PATH, LIBRARIES_TO_BE_PARSED_EXTEND,
                  OPENMS_GIT_LC_DATE_FORMAT, OPENMP_FOUND, OPENMP_CXX_FLAGS)

IS_DEBUG = OPEN_MS_BUILD_TYPE.upper() == "DEBUG"
OMP = (OPENMP_FOUND.upper() == "ON" or OPENMP_FOUND.upper() == "TRUE" or OPENMP_FOUND == "1")

if iswin and IS_DEBUG:
    raise Exception("building pyopenms on windows in debug mode not tested yet.")

# use autowrap to generate Cython and .cpp file for wrapping OpenMS:
import pickle
import os
import glob
import re
import shutil
import time

if OPEN_MS_GIT_BRANCH == "nightly":
    package_name = "pyopenms"
    package_version = OPEN_MS_VERSION + ".dev" + OPENMS_GIT_LC_DATE_FORMAT
else:
    package_name = "pyopenms"
    package_version = OPEN_MS_VERSION

os.environ["CC"] = OPEN_MS_COMPILER
# AFAIK distutils does not care about CXX (set it to be safe)
os.environ["CXX"] = OPEN_MS_COMPILER

j = os.path.join

src_pyopenms = PYOPENMS_SRC_DIR
extra_includes = glob.glob(src_pyopenms + "/extra_includes/*.h*")

for include in extra_includes:
    shutil.copy(include, "extra_includes/")


persisted_data_path = "include_dir.bin"
autowrap_include_dirs = pickle.load(open(persisted_data_path, "rb"))

from setuptools import setup, Extension

with open("pyopenms/version.py", "w") as fp:
    print("version=%r" % package_version, file=fp)

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
    if IS_DEBUG:
        libraries = ["OpenMSd", "OpenSwathAlgod", "Qt5Cored", "Qt5Networkd"]
    else:
        libraries = ["OpenMS", "OpenSwathAlgo", "Qt5Core", "Qt5Network"]
elif sys.platform.startswith("linux"):
    libraries = ["OpenMS", "OpenSwathAlgo", "Qt5Core", "Qt5Network"]
elif sys.platform == "darwin":
    libraries = ["OpenMS", "OpenSwathAlgo"]
else:
    print("\n")
    print("platform", sys.platform, "not supported yet")
    print("\n")
    exit()

if (iswin):
    library_dirs = [OPEN_MS_BUILD_DIR,
                    j(OPEN_MS_BUILD_DIR, "lib", "Release"),
                    j(OPEN_MS_BUILD_DIR, "bin", "Release"),
                    j(OPEN_MS_BUILD_DIR, "Release"),
                    QT_INSTALL_BINS,
                    QT_INSTALL_LIBS,
                    ]
else:
    library_dirs = [OPEN_MS_BUILD_DIR,
                j(OPEN_MS_BUILD_DIR, "lib"),
                j(OPEN_MS_BUILD_DIR, "bin"),
                QT_INSTALL_BINS,
                QT_INSTALL_LIBS,
                ]

# extend with contrib lib dirs
for OPEN_MS_CONTRIB_BUILD_DIR in OPEN_MS_CONTRIB_BUILD_DIRS.split(";"):
  library_dirs.append(j(OPEN_MS_CONTRIB_BUILD_DIR, "lib"))

import numpy


include_dirs = [
    "extra_includes",
    j(numpy.core.__path__[0], "include")
]

# append all include and library dirs exported by CMake
include_dirs.extend(PYOPENMS_INCLUDE_DIRS.split(";"))

if INCLUDE_DIRS_EXTEND: # only add if not empty
    include_dirs.extend(INCLUDE_DIRS_EXTEND.split(";"))
if LIBRARY_DIRS_EXTEND: # only add if not empty
    library_dirs.extend(LIBRARY_DIRS_EXTEND.split(";"))
if LIBRARIES_EXTEND: # only add if not empty
    libraries.extend(LIBRARIES_EXTEND.split(";"))


# libraries of any type to be parsed and added
objects = []
add_libs = LIBRARIES_TO_BE_PARSED_EXTEND.split(";")
for lib in add_libs:
  if not iswin:
    if lib.endswith(".a"):
      objects.append(lib)
      name_search = re.search('.*/lib(.*)\.a$', lib)
      if name_search:
        libraries.append(name_search.group(1))
        library_dirs.append(os.path.dirname(lib))
    if lib.endswith(".so") or lib.endswith(".dylib"):
      name_search = re.search('.*/lib(.*)\.(so|dylib)$', lib)
      if name_search:
        libraries.append(name_search.group(1))
        library_dirs.append(os.path.dirname(lib))
  else:
    if lib.endswith(".lib"):
      name_search = re.search('.*/(.*)\.lib$', lib)
      if name_search:
        libraries.append(name_search.group(1))
        library_dirs.append(os.path.dirname(lib))


extra_link_args = []
extra_compile_args = []

if iswin:
    # /EHs is important. It sets _CPPUNWIND which causes boost to
    # set BOOST_NO_EXCEPTION in <boost/config/compiler/visualc.hpp>
    # such that  boost::throw_excption() is declared but not implemented.
    # The linker does not like that very much ...
    extra_compile_args = ["/EHs", "/bigobj"]
    extra_compile_args.append("/std:c++17")
    extra_link_args.append("/std:c++17")

elif sys.platform.startswith("linux"):
    extra_link_args = ["-Wl,-s"]
    if OMP:
        libraries.append("gomp")
        libraries.append("pthread")
elif sys.platform == "darwin":
    library_dirs.insert(0,j(OPEN_MS_BUILD_DIR,"pyOpenMS","pyopenms"))
    if OMP:
        libraries.append("omp")
    # we need to manually link to the Qt Frameworks
    extra_compile_args = ["-Qunused-arguments"]
    extra_link_args = ["-Wl,-rpath","-Wl,@loader_path/"]
if IS_DEBUG:
    extra_compile_args.append("-g2")
if OMP and OPENMP_CXX_FLAGS:
    extra_compile_args.extend(OPENMP_CXX_FLAGS.split(";"))

if not iswin:
    extra_link_args.append("-std=c++17")
    extra_compile_args.append("-std=c++17")
    if isosx: # MacOS
        extra_compile_args.append("-stdlib=libc++")
        extra_link_args.append("-stdlib=libc++") # MacOS libstdc++ does not include c++11+ lib support.
        extra_link_args.append("-mmacosx-version-min=10.9") # due to libc++
        extra_compile_args.append("-Wno-deprecated")
        extra_compile_args.append("-Wno-nullability-completeness")
        if (osx_ver >= "10.14.0" and SYSROOT_OSX_PATH): # since macOS Mojave
            extra_link_args.append("-isysroot" + SYSROOT_OSX_PATH)
            extra_compile_args.append("-isysroot" + SYSROOT_OSX_PATH)
    else:
        extra_compile_args.append("-Wno-deprecated-copy")
    extra_compile_args.append("-Wno-redeclared-class-member")
    extra_compile_args.append("-Wno-unused-local-typedefs")
    extra_compile_args.append("-Wdeprecated-declarations")
    extra_compile_args.append("-Wno-sign-compare")
    extra_compile_args.append("-Wno-unknown-pragmas")
    extra_compile_args.append("-Wno-header-guard")
    extra_compile_args.append("-Wno-unused-function")
    extra_compile_args.append("-Wno-deprecated-declarations")
    extra_compile_args.append("-Wno-missing-declarations")
    extra_compile_args.append("-Wno-int-in-bool-context")
    
    if no_optimization:
        extra_compile_args.append("-O0")
        extra_link_args.append("-O0")

mnames = ["_pyopenms_%s" % (k+1) for k in range(int(PY_NUM_MODULES))]
ext = []

##WARNING debug
libraries.extend("boost_regex-mt-x64")

for module in mnames:

    ext.append(Extension(
        module,
        sources=["pyopenms/%s.cpp" % module],
        language="c++",
        library_dirs=library_dirs,
        libraries=libraries,
        include_dirs=include_dirs + autowrap_include_dirs,
        extra_compile_args=extra_compile_args,
        extra_objects=objects,
        extra_link_args=extra_link_args,
		define_macros=[('BOOST_ALL_NO_LIB', None), ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")] ## Deactivates boost autolink (esp. on win). Shuts up the damn "deprecated NumPy API" warning spam (https://docs.cython.org/en/latest/src/userguide/numpy_tutorial.html#numpy-compilation)
		## Alternative is to specify the boost naming scheme (--layout param; easy if built from contrib)
		## TODO just take over compile definitions from OpenMS (CMake)
    ))

# enforce 64bit-only build as OpenMS is not available in 32bit on osx
if sys.platform == "darwin":
    os.environ['ARCHFLAGS'] = "-arch x86_64"

setup(

    name=package_name,
    packages=["pyopenms"],
    ext_package="pyopenms",
    package_data= {
        'pyopenms': ['py.typed', '*.pyi']
    },
	install_requires=[
          'numpy<=1.26.4',
          'pandas',
          'matplotlib>=3.5'
    ],

    version=package_version,

    maintainer="The OpenMS team",
    maintainer_email="open-ms-general@lists.sourceforge.net",
    license="http://opensource.org/licenses/BSD-3-Clause",
    platforms=["any"],
    description="Python wrapper for C++ LC-MS library OpenMS",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: BSD License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    long_description=open("README.rst").read(),
    long_description_content_type="text/x-rst",
    zip_safe=False,

    url="https://openms.de",
    project_urls={
        "Documentation": "https://pyopenms.readthedocs.io",
        "Source Code": "https://github.com/OpenMS/OpenMS/tree/develop/src/pyOpenMS",
        "Tracker": "https://github.com/OpenMS/OpenMS/issues",
        "Documentation Source": "https://github.com/OpenMS/pyopenms-docs",
    },

    author="OpenMS team",
    author_email="webmaster@openms.de",

    ext_modules=ext,
    include_package_data=True  # see MANIFEST.in
)
