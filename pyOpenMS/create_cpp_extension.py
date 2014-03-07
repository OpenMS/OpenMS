# input-encoding: latin-1

# windows ?
import sys
iswin = sys.platform == "win32"

# import config
from env import (QT_QMAKE_VERSION_INFO, QT_LIBRARY_DIR, OPEN_MS_BUILD_TYPE, OPEN_MS_SRC,
                 OPEN_MS_CONTRIB_BUILD_DIRS, OPEN_MS_LIB, OPEN_SWATH_ALGO_LIB, OPEN_MS_BUILD_DIR,
                 MSVCR90DLL, MSVCP90DLL, OPEN_MS_VERSION)

IS_DEBUG = OPEN_MS_BUILD_TYPE.upper() == "DEBUG"

if iswin and IS_DEBUG:
    raise Exception("building pyopenms on windows in debug mode not tested yet.")

# use autowrap to generate Cython and .cpp file for wrapping OpenMS:
import autowrap.Main
import glob
import cPickle
import os.path
import os
import shutil
import time

j = os.path.join

src_pyopenms = j(OPEN_MS_SRC, "pyOpenMS")
pxd_files = glob.glob(src_pyopenms + "/pxds/*.pxd")
addons = glob.glob(src_pyopenms + "/addons/*.pyx")
converters = [j(src_pyopenms, "converters")]


persisted_data_path = "include_dir.bin"

extra_cimports = [  # "from libc.stdint cimport *",
    #"from libc.stddef cimport *",
    #"from UniqueIdInterface cimport setUniqueId as _setUniqueId",
    #"from Map cimport Map as _Map",
    #"cimport numpy as np"
]

decls, instance_map = autowrap.parse(pxd_files, ".")
# add __str__ if toString() method is declared:
for d in decls:
    # enums, free functions, .. do not have a methods attribute
    methods = getattr(d, "methods", dict())
    to_strings = []
    for name, mdecls in methods.items():
        for mdecl in mdecls:
            name = mdecl.cpp_decl.annotations.get("wrap-cast", name)
            name = mdecl.cpp_decl.annotations.get("wrap-as", name)
            if name == "toString":
                to_strings.append(mdecl)

    for to_string in to_strings:
        if len(to_string.arguments) == 0:
            d.methods.setdefault("__str__", []).append(to_string)
            print "ADDED __str__ method to", d.name
            break

autowrap_include_dirs = autowrap.Main.create_wrapper_code(decls, instance_map, addons,
                                                          converters, "pyopenms/pyopenms.pyx",
                                                          extra_cimports,
                                                          None)

cPickle.dump(autowrap_include_dirs, open(persisted_data_path, "wb"))

import re
f = open("pyopenms/pyopenms.cpp")
fout = open("pyopenms/pyopenms_out.cpp", "w")
expr_fix = re.compile(r"(.*).std::vector<(.*)>::iterator::~iterator\(\)")
for line in f:
    res = expr_fix.sub('typedef std::vector<\\2>::iterator _it;\n\\1.~_it()', line)
    fout.write(res)

fout.close()
f.close()
shutil.copy("pyopenms/pyopenms_out.cpp", "pyopenms/pyopenms.cpp")
os.remove("pyopenms/pyopenms_out.cpp")

print "created pyopenms.cpp"

# create version information
version = OPEN_MS_VERSION

print >> open("pyopenms/version.py", "w"), "version=%r\n" % version
print >> open("pyopenms/qt_version_info.py", "w"), "info=%r\n" % QT_QMAKE_VERSION_INFO

# parse config

if OPEN_MS_CONTRIB_BUILD_DIRS.endswith(";"):
    OPEN_MS_CONTRIB_BUILD_DIRS = OPEN_MS_CONTRIB_BUILD_DIRS[:-1]

for OPEN_MS_CONTRIB_BUILD_DIR in OPEN_MS_CONTRIB_BUILD_DIRS.split(";"):
    if os.path.exists(os.path.join(OPEN_MS_CONTRIB_BUILD_DIR, "lib")):
        break


if iswin:
    # fix for broken library names in release 1.11:
    for p in glob.glob(os.path.join(OPEN_MS_CONTRIB_BUILD_DIR,
                                    "lib",
                                    "libboost_math_*mt.lib")):

        if "vc90" in p:
            continue
        new_p = p.replace("-mt.lib", "-vc90-mt-1_52.lib")
        shutil.copy(p, new_p)


# Package data expected to be installed. On Linux the debian package
# contains share/ data and must be installed to get access to the OpenMS shared
# library.
#
if iswin:
    shutil.copy(OPEN_MS_LIB, "pyopenms")
    shutil.copy(OPEN_SWATH_ALGO_LIB, "pyopenms")

    shutil.copy(MSVCR90DLL, "pyopenms")
    shutil.copy(MSVCP90DLL, "pyopenms")

    if OPEN_MS_BUILD_TYPE.upper() == "DEBUG":
        shutil.copy(j(QT_LIBRARY_DIR, "QtCored4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtGuid4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtSqld4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtNetworkd4.dll"), "pyopenms")
        shutil.copy(j(OPEN_MS_CONTRIB_BUILD_DIR, "lib", "xerces-c_3_1D.dll"),
                    "pyopenms")
    else:
        shutil.copy(j(QT_LIBRARY_DIR, "QtCore4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtGui4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtSql4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtNetwork4.dll"), "pyopenms")
        shutil.copy(j(OPEN_MS_CONTRIB_BUILD_DIR, "lib", "xerces-c_3_1.dll"),
                    "pyopenms")

elif sys.platform == "linux2":

    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenMS.so"), "pyopenms")
    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenSwathAlgo.so"), "pyopenms")

elif sys.platform == "darwin":

    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenMS.dylib"), "pyopenms")
    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenSwathAlgo.dylib"), "pyopenms")

else:
    print
    print "platform ", sys.platform, "not supported yet"
    print
    exit()

print "copied files needed for distribution to pyopenms/"
print

