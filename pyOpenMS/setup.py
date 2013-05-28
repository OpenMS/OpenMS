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
import os, shutil

j=os.path.join

pxd_files = glob.glob("pxds/*.pxd")
addons = glob.glob("addons/*.pyx")
converters = ["converters"]

def copy_files_if_updated(src, target, ext):
    getf = lambda d: [f for f in os.listdir(d) if os.path.splitext(f)[1]==ext]
    files_target = getf(src)
    files_src = getf(target)

    for f in files_src:
        src_f = j(src, f)
        if f not in files_target:
            print "COPY NEW", src_f
            shutil.copy(src_f, target)
        else:
            src_time = os.path.getmtime(src_f)
            target_time = os.path.getmtime(j(target, f))
            if src_time > target_time:
                print "UPDATE", src_f
                shutil.copy(src_f, target)


pp = j(OPEN_MS_SRC, "pyOpenMS")
copy_files_if_updated(j(pp, "pxds"), "pxds", ".pxd")
copy_files_if_updated(j(pp, "addons"), "addons", ".pyx")
copy_files_if_updated(j(pp, "converters"), "converters", ".py")

# sometimes I use symlinks for restricting pxd folder and others, so
# resolve, else os.path.getmtime will not work
real_pathes = [os.path.realpath(p) for p in pxd_files + addons + converters]
mtimes = [os.path.getmtime(f) for f in real_pathes]

if os.path.exists("pyopenms/pyopenms.pyx"):
    mtime_result = os.path.getmtime("pyopenms/pyopenms.pyx")
else:
    mtime_result = 0

persisted_data_path = "include_dir.bin"

if not os.path.exists(persisted_data_path)\
    or any(m > mtime_result for m in mtimes):

    extra_cimports = [ # "from libc.stdint cimport *",
                #"from libc.stddef cimport *",
                #"from UniqueIdInterface cimport setUniqueId as _setUniqueId",
                "from Map cimport Map as _Map",
                #"cimport numpy as np"
                ]
    autowrap_include_dirs = autowrap.Main.run(pxd_files,
                                            addons,
                                            converters,
                                            "pyopenms/pyopenms.pyx",
                                            extra_cimports)
    cPickle.dump(autowrap_include_dirs, open(persisted_data_path, "wb"))

else:
    autowrap_include_dirs = cPickle.load(open(persisted_data_path, "rb"))

from setuptools import setup, Extension
import time

# Due to a bug in Cython when dealing with destructors of typedefs, we have to
# fix the .cpp file manually (the bug is only triggered with clang, see also
# https://bugzilla.mozilla.org/show_bug.cgi?id=623303 )
# We thus replace all occurences of
#   p->__pyx_v_it.std::vector<T>::iterator::~iterator();
# with
#   typedef std::vector<T>::iterator _it;
#   p->__pyx_v_it.~_it();

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

# create version information

ctime = os.stat("pyopenms").st_mtime
ts = time.gmtime(ctime)
timestamp = "%02d-%02d-%4d" % (ts.tm_mday, ts.tm_mon, ts.tm_year)

from version import version

print >> open("pyopenms/version.py", "w"), "version=%r\n" % version
print >> open("pyopenms/qt_version_info.py", "w"), "info=%r\n" % QT_QMAKE_VERSION_INFO


# parse config

if OPEN_MS_CONTRIB_BUILD_DIRS.endswith(";"):
    OPEN_MS_CONTRIB_BUILD_DIRS = OPEN_MS_CONTRIB_BUILD_DIRS[:-1]

for OPEN_MS_CONTRIB_BUILD_DIR in  OPEN_MS_CONTRIB_BUILD_DIRS.split(";"):
    if os.path.exists(OPEN_MS_CONTRIB_BUILD_DIR):
        break



# package data expected to be installed. on linux the debian package
# contains share/ data and must be intalled to get access to the openms shared
# library.
#
if iswin:
    shutil.copy(OPEN_MS_LIB, "pyopenms")
    shutil.copy(OPEN_SWATH_ALGO_LIB, "pyopenms")

    shutil.copy(MSVCR90DLL, "pyopenms")
    shutil.copy(MSVCP90DLL, "pyopenms")


    if OPEN_MS_BUILD_TYPE.upper() == "DEBUG":
        libraries=["OpenMSd", "OpenSwathAlgod", "xerces-c_3D", "QtCored4", "gsl_d", "cblas_d"]
        shutil.copy(j(QT_LIBRARY_DIR, "QtCored4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtGuid4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtSqld4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtNetworkd4.dll"), "pyopenms")
        shutil.copy(j(OPEN_MS_CONTRIB_BUILD_DIR, "lib", "xerces-c_3_1D.dll"),\
                        "pyopenms")
    else:
        libraries=["OpenMS", "OpenSwathAlgo", "xerces-c_3", "QtCore4", "gsl", "cblas"]
        shutil.copy(j(QT_LIBRARY_DIR, "QtCore4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtGui4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtSql4.dll"), "pyopenms")
        shutil.copy(j(QT_LIBRARY_DIR, "QtNetwork4.dll"), "pyopenms")
        shutil.copy(j(OPEN_MS_CONTRIB_BUILD_DIR, "lib", "xerces-c_3_1.dll"),\
                        "pyopenms")

elif sys.platform == "linux2":
    libraries=["OpenMS", "OpenSwathAlgo", "xerces-c", "QtCore", "gsl", "gslcblas"]

    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenMS.so" ), "pyopenms")
    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenSwathAlgo.so" ), "pyopenms")

elif sys.platform == "darwin":
    libraries=["OpenMS", "OpenSwathAlgo", "xerces-c", "gsl", "gslcblas" ]

    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenMS.dylib" ), "pyopenms")
    shutil.copy(j(OPEN_MS_BUILD_DIR, "lib", "libOpenSwathAlgo.dylib" ), "pyopenms")

else:
    print
    print "platform ", sys.platform, "not supported yet"
    print
    exit()

library_dirs=[ OPEN_MS_BUILD_DIR,
               j(OPEN_MS_BUILD_DIR,"lib"),
               j(OPEN_MS_BUILD_DIR,"bin"),
               j(OPEN_MS_BUILD_DIR,"bin", "Release"),
               j(OPEN_MS_BUILD_DIR,"Release"),
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

extra_link_args = []

if sys.platform == "linux2":
    extra_link_args = [ "-Wl,-s"]
elif sys.platform == "darwin":
    # we need to manually link to the Qt Frameworks
    extra_link_args = [ "-Wl,-s",
                        "-F"+QT_LIBRARY_DIR,
                        "-framework Carbon",
                        "-framework AGL",
                        "-framework OpenGL",
                        "-framework QtOpenGL",
                        "-framework QtSvg",
                        "-framework QtWebKit",
                        "-framework QtXmlPatterns",
                        "-framework QtGui",
                        "-framework QtTest",
                        "-framework QtXml",
                        "-framework QtSql",
                        "-framework QtNetwork",
                        "-framework QtCore" ]

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
        extra_compile_args = iswin and [ "/EHs", "/bigobj"] or (IS_DEBUG and ["-g2"] or []),
        extra_link_args = extra_link_args
    )


source_share = j(OPEN_MS_SRC, "share")
target_share = "pyopenms/share"
if os.path.exists(target_share):
    shutil.rmtree(target_share)
shutil.copytree(source_share, target_share, ignore=lambda *a: ["examples"])

share_data = []

if iswin:
    share_data += [MSVCR90DLL, "xerces-c_3_1.dll"]


share_data.append("License.txt")

# enforce 64bit-only build as OpenMS is not available in 32bit on osx
if sys.platform == "darwin":
    os.environ['ARCHFLAGS'] = "-arch x86_64"

setup(

    name = "pyopenms",
    packages = ["pyopenms"],
    ext_package = "pyopenms",

    version = version,

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

    url = "http://open-ms.de",

    author = "Uwe Schmitt",
    author_email = "uschmitt@mineway.de",

    ext_modules = [ext ],
    include_package_data = True # see MANIFEST.in
)
